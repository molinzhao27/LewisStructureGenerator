#include "lewis_engine.h"

#include <string.h>

typedef struct {
    uint8_t valence_pairs;
    uint8_t bond_pairs;
    uint8_t lone_pairs;
    const char *ep_geometry;
    const char *shape;
    const char *hybridization;
    const char *bond_angle;
} VseprRow;

/* Rows with a single bonded atom describe diatomics, which have no bond angle. */
static const VseprRow vsepr_rows[] = {
    {1, 1, 0, "Linear", "Linear", "s", "N/A"},
    {2, 2, 0, "Linear", "Linear", "sp", "180"},
    {2, 1, 1, "Linear", "Linear", "sp", "N/A"},
    {3, 3, 0, "Trigonal Planar", "Trigonal Planar", "sp2", "120"},
    {3, 2, 1, "Trigonal Planar", "Bent", "sp2", "<120"},
    {3, 1, 2, "Trigonal Planar", "Linear", "sp2", "N/A"},
    {4, 4, 0, "Tetrahedral", "Tetrahedral", "sp3", "109.5"},
    {4, 3, 1, "Tetrahedral", "Trigonal Pyramidal", "sp3", "<109.5"},
    {4, 2, 2, "Tetrahedral", "Bent", "sp3", "<109.5"},
    {4, 1, 3, "Tetrahedral", "Linear", "sp3", "N/A"},
    {5, 5, 0, "Trigonal Bipyramidal", "Trigonal Bipyramidal", "sp3d", "90, 120"},
    {5, 4, 1, "Trigonal Bipyramidal", "Seesaw", "sp3d", "<90, <120"},
    {5, 3, 2, "Trigonal Bipyramidal", "T-shaped", "sp3d", "<90"},
    {5, 2, 3, "Trigonal Bipyramidal", "Linear", "sp3d", "180"},
    {5, 1, 4, "Trigonal Bipyramidal", "Linear", "sp3d", "N/A"},
    {6, 6, 0, "Octahedral", "Octahedral", "sp3d2", "90"},
    {6, 5, 1, "Octahedral", "Square Pyramidal", "sp3d2", "<90"},
    {6, 4, 2, "Octahedral", "Square Planar", "sp3d2", "90"},
    {6, 3, 3, "Octahedral", "T-shaped", "sp3d2", "<90"},
    {6, 2, 4, "Octahedral", "Linear", "sp3d2", "180"},
    {6, 1, 5, "Octahedral", "Linear", "sp3d2", "N/A"},
    {7, 7, 0, "Pentagonal Bipyramidal", "Pentagonal Bipyramidal", "sp3d3", "72, 90"},
    {7, 6, 1, "Pentagonal Bipyramidal", "Pentagonal Pyramidal", "sp3d3", "<90"},
    {7, 5, 2, "Pentagonal Bipyramidal", "Pentagonal Planar", "sp3d3", "72"},
    {7, 4, 3, "Pentagonal Bipyramidal", "Seesaw", "sp3d3", "<90, <72"},
    {7, 3, 4, "Pentagonal Bipyramidal", "T-shaped", "sp3d3", "<90"},
    {7, 2, 5, "Pentagonal Bipyramidal", "Linear", "sp3d3", "180"},
    {7, 1, 6, "Pentagonal Bipyramidal", "Linear", "sp3d3", "N/A"},
};

/*
 * Structure generation
 *
 * Every plausible sigma skeleton (a tree) is built, and for each one every
 * distribution of pi bonds over its heavy-atom bonds is evaluated. Lone pairs
 * complete octets on the most electronegative atoms first; any electrons left
 * over go to an expandable (period 3+) atom. Forms are ranked octet-first:
 *
 *   1. fewest electrons beyond an octet (expand only when unavoidable)
 *   2. fewest electrons missing from carbon-group octets (carbocations)
 *   3. smallest total |formal charge|
 *
 * Every form of the winning skeleton that ties on those keys is a resonance
 * structure. Forms are ordered to avoid concentrated charges (N=N=N before
 * N#N-N for azide), then so negative charge sits on the more electronegative
 * atoms. Between skeletons that tie, more resonance forms wins (allyl cation).
 *
 * Search state is file-static to keep stack use small on the calculator.
 */

typedef struct {
    int excess;       /* electrons beyond an octet */
    int deficit;      /* electrons missing from group-14 octets */
    int sum_abs_fc;   /* sum of |formal charge| */
    int sum_sq_fc;    /* sum of formal charge squared */
    int en_charge;    /* sum of formal charge * electronegativity */
} FormScore;

typedef struct {
    FormScore form;          /* score of the skeleton's best form */
    uint8_t num_forms;
    int max_heavy_degree;    /* prefer skeletons that cluster atoms on one center */
    int h_host_en;           /* prefer H on less electronegative atoms */
} SkeletonScore;

static Molecule *g_mol;
static uint8_t g_n;
static uint8_t g_elem[MAX_ATOMS];
static int g_ve;

static uint8_t g_heavy[MAX_ATOMS];
static uint8_t g_n_heavy;
static uint8_t g_hydrogen[MAX_ATOMS];
static uint8_t g_n_hydrogen;
static uint8_t g_fill_order[MAX_ATOMS];   /* atoms by decreasing electronegativity */

/* Skeleton under construction */
static LewisStructure g_ls;
static uint8_t g_sigma[MAX_ATOMS];
static uint8_t g_bos[MAX_ATOMS];
static uint8_t g_hh_bonds[MAX_BONDS];
static uint8_t g_n_hh;
static int8_t g_expand_atom;

static uint8_t g_backbone[MAX_ATOMS];
static uint8_t g_n_backbone;
static uint8_t g_pool[MAX_ATOMS];
static bool g_pool_used[MAX_ATOMS];
static uint8_t g_n_pool;
static uint8_t g_terminal[MAX_ATOMS];
static uint8_t g_n_terminal;
static uint8_t g_term_host[MAX_ATOMS];
static uint8_t g_h_cap[MAX_ATOMS];
static uint8_t g_h_count[MAX_ATOMS];

/* Forms of the current skeleton, best first */
static LewisStructure g_cur[MAX_RESONANCE];
static FormScore g_cur_score[MAX_RESONANCE];
static uint8_t g_cur_n;

/* Best skeleton so far; its forms live in g_mol->res */
static SkeletonScore g_best_score;
static uint8_t g_best_n;

static InvalidReason g_fail;

static int abs_int(int x)
{
    return (x < 0) ? -x : x;
}

static const Element *elem_of(uint8_t atom)
{
    return &elements[g_elem[atom]];
}

static bool is_duet(uint8_t atom)
{
    return g_elem[atom] == ELEM_H || g_elem[atom] == ELEM_HE;
}

/* Groups 2 and 13 (Be, B, Al...) may keep incomplete octets and hold no lone pairs. */
static bool is_sigma_only(uint8_t atom)
{
    uint8_t group = elem_of(atom)->group;
    return group == 2 || group == 13;
}

/* Bonds formed by the uncharged atom: C 4, N 3, O 2, F 1. */
static int natural_valence(uint8_t elem_idx)
{
    const Element *e = &elements[elem_idx];
    if (e->group == 18) return 0;
    return (e->valence <= 4) ? e->valence : 8 - e->valence;
}

static uint8_t max_sigma(uint8_t elem_idx)
{
    const Element *e = &elements[elem_idx];
    if (elem_idx == ELEM_H) return 1;
    if (e->bond_cap == 0) return 0;          /* He, Ne, Ar */
    if (e->group == 1) return 1;
    if (e->period == 2) {
        if (e->group == 17) return 1;
        if (e->group == 16) return 3;        /* H3O+ */
        return 4;                            /* NH4+, BF4- */
    }
    if (e->group == 2) return 4;
    return (e->group >= 17) ? 7 : 6;
}

static uint8_t max_bond_order_sum(uint8_t atom)
{
    if (g_elem[atom] == ELEM_H) return 1;
    return (elem_of(atom)->period <= 2) ? 4 : 7;
}

/* Atoms that chain together to form a backbone (B, C, N, Si, P...). */
static bool is_backbone_elem(uint8_t elem_idx)
{
    uint8_t group = elements[elem_idx].group;
    return group >= 13 && group <= 15;
}

static bool can_expand(uint8_t atom)
{
    const Element *e = elem_of(atom);
    return g_elem[atom] != ELEM_H && e->period >= 3 && e->group >= 14;
}

static void note_failure(InvalidReason reason)
{
    /* Keep the reason from the furthest stage reached. */
    static const uint8_t rank[] = {
        0, /* INVALID_NONE */
        0, /* INVALID_NO_ATOMS */
        0, /* INVALID_NEGATIVE_ELECTRONS */
        0, /* INVALID_ODD_ELECTRONS */
        1, /* INVALID_SKELETON */
        2, /* INVALID_LEFTOVER_ELECTRONS */
        3, /* INVALID_SHELL_RULE */
        0  /* INVALID_FORMAL_CHARGE_SUM */
    };
    if (rank[reason] > rank[g_fail]) {
        g_fail = reason;
    }
}

static bool push_bond(uint8_t a, uint8_t b)
{
    if (g_ls.num_bonds >= MAX_BONDS) return false;
    if (g_sigma[a] >= max_sigma(g_elem[a])) return false;
    if (g_sigma[b] >= max_sigma(g_elem[b])) return false;

    Bond *bd = &g_ls.bonds[g_ls.num_bonds++];
    bd->a = a;
    bd->b = b;
    bd->order = 1;
    g_sigma[a]++;
    g_sigma[b]++;
    return true;
}

static void pop_bond(void)
{
    Bond *bd = &g_ls.bonds[--g_ls.num_bonds];
    g_sigma[bd->a]--;
    g_sigma[bd->b]--;
}

/* Negative: a ranks better. Compares the octet-first keys only. */
static int compare_rank(const FormScore *a, const FormScore *b)
{
    if (a->excess != b->excess) return (a->excess < b->excess) ? -1 : 1;
    if (a->deficit != b->deficit) return (a->deficit < b->deficit) ? -1 : 1;
    if (a->sum_abs_fc != b->sum_abs_fc) return (a->sum_abs_fc < b->sum_abs_fc) ? -1 : 1;
    return 0;
}

/* Negative: a should be listed first among equally ranked forms. */
static int compare_order(const FormScore *a, const FormScore *b)
{
    if (a->sum_sq_fc != b->sum_sq_fc) return (a->sum_sq_fc < b->sum_sq_fc) ? -1 : 1;
    if (a->en_charge != b->en_charge) return (a->en_charge < b->en_charge) ? -1 : 1;
    return 0;
}

static void record_form(const FormScore *s)
{
    if (g_cur_n > 0) {
        int c = compare_rank(s, &g_cur_score[0]);
        if (c > 0) return;
        if (c < 0) g_cur_n = 0;
    }

    uint8_t pos = g_cur_n;
    while (pos > 0 && compare_order(&g_cur_score[pos - 1], s) > 0) {
        pos--;
    }
    if (pos >= MAX_RESONANCE) return;

    uint8_t last = (g_cur_n < MAX_RESONANCE) ? g_cur_n : (MAX_RESONANCE - 1);
    for (uint8_t k = last; k > pos; k--) {
        g_cur[k] = g_cur[k - 1];
        g_cur_score[k] = g_cur_score[k - 1];
    }
    g_cur[pos] = g_ls;
    g_cur_score[pos] = *s;
    if (g_cur_n < MAX_RESONANCE) g_cur_n++;
}

/* Place lone pairs for the current bond orders, then score the form. */
static void evaluate_form(void)
{
    int pairs = g_ve / 2;
    for (uint8_t b = 0; b < g_ls.num_bonds; b++) {
        pairs -= g_ls.bonds[b].order;
    }
    if (pairs < 0) {
        note_failure(INVALID_SHELL_RULE);
        return;
    }

    memset(g_ls.lone_pairs, 0, sizeof(g_ls.lone_pairs));
    for (uint8_t k = 0; k < g_n && pairs > 0; k++) {
        uint8_t i = g_fill_order[k];
        if (is_duet(i) || is_sigma_only(i)) continue;
        int need = 4 - g_bos[i];
        if (need <= 0) continue;
        if (need > pairs) need = pairs;
        g_ls.lone_pairs[i] = (uint8_t)need;
        pairs -= need;
    }
    if (pairs > 0) {
        if (g_expand_atom < 0) {
            note_failure(INVALID_LEFTOVER_ELECTRONS);
            return;
        }
        g_ls.lone_pairs[g_expand_atom] = (uint8_t)(g_ls.lone_pairs[g_expand_atom] + pairs);
    }

    FormScore s = {0, 0, 0, 0, 0};
    for (uint8_t i = 0; i < g_n; i++) {
        const Element *e = elem_of(i);
        int electrons = 2 * (g_ls.lone_pairs[i] + g_bos[i]);

        if (is_duet(i)) {
            if (electrons != 2) {
                note_failure(INVALID_SHELL_RULE);
                return;
            }
        } else {
            if (electrons > ((e->period <= 2) ? 8 : 14)) {
                note_failure(INVALID_SHELL_RULE);
                return;
            }
            if (electrons > 8) {
                s.excess += electrons - 8;
            } else if (electrons < 8 && !is_sigma_only(i)) {
                if (e->group != 14) {
                    note_failure(INVALID_SHELL_RULE);
                    return;
                }
                s.deficit += 8 - electrons;
            }
        }

        int fc = e->valence - 2 * g_ls.lone_pairs[i] - g_bos[i];
        g_ls.formal_charge[i] = (int8_t)fc;
        s.sum_abs_fc += abs_int(fc);
        s.sum_sq_fc += fc * fc;
        s.en_charge += fc * e->eneg;
    }

    record_form(&s);
}

static void enumerate_pi(uint8_t k, int remaining)
{
    if (remaining > 2 * (g_n_hh - k)) return;
    if (k == g_n_hh) {
        evaluate_form();
        return;
    }

    Bond *bd = &g_ls.bonds[g_hh_bonds[k]];
    int max_extra = (remaining < 2) ? remaining : 2;
    for (int extra = max_extra; extra >= 0; extra--) {
        if (g_bos[bd->a] + extra > max_bond_order_sum(bd->a)) continue;
        if (g_bos[bd->b] + extra > max_bond_order_sum(bd->b)) continue;

        bd->order = (uint8_t)(1 + extra);
        g_bos[bd->a] = (uint8_t)(g_bos[bd->a] + extra);
        g_bos[bd->b] = (uint8_t)(g_bos[bd->b] + extra);
        enumerate_pi(k + 1, remaining - extra);
        g_bos[bd->a] = (uint8_t)(g_bos[bd->a] - extra);
        g_bos[bd->b] = (uint8_t)(g_bos[bd->b] - extra);
        bd->order = 1;
    }
}

static bool skeleton_better(const SkeletonScore *a, const SkeletonScore *b)
{
    int c = compare_rank(&a->form, &b->form);
    if (c != 0) return c < 0;
    if (a->num_forms != b->num_forms) return a->num_forms > b->num_forms;
    if (a->max_heavy_degree != b->max_heavy_degree) return a->max_heavy_degree > b->max_heavy_degree;
    if (a->h_host_en != b->h_host_en) return a->h_host_en < b->h_host_en;
    return compare_order(&a->form, &b->form) < 0;
}

static void evaluate_skeleton(void)
{
    uint8_t n_sigma = g_ls.num_bonds;
    if (g_ve < 2 * n_sigma) {
        note_failure(INVALID_SKELETON);
        return;
    }

    /* Bond order total that gives every atom exactly its target shell. */
    int target = 0;
    for (uint8_t i = 0; i < g_n; i++) {
        g_bos[i] = g_sigma[i];
        if (is_duet(i)) target += 2;
        else if (is_sigma_only(i)) target += 2 * g_sigma[i];
        else target += 8;
    }
    int max_pi = (target - g_ve) / 2 - n_sigma;
    if (max_pi < 0) max_pi = 0;

    g_n_hh = 0;
    g_expand_atom = -1;
    for (uint8_t b = 0; b < n_sigma; b++) {
        if (g_elem[g_ls.bonds[b].a] != ELEM_H && g_elem[g_ls.bonds[b].b] != ELEM_H) {
            g_hh_bonds[g_n_hh++] = b;
        }
    }
    for (uint8_t i = 0; i < g_n; i++) {
        if (!can_expand(i)) continue;
        if (g_expand_atom < 0 ||
            g_sigma[i] > g_sigma[g_expand_atom] ||
            (g_sigma[i] == g_sigma[g_expand_atom] && elem_of(i)->eneg < elem_of(g_expand_atom)->eneg)) {
            g_expand_atom = (int8_t)i;
        }
    }

    /* More pi bonds than max_pi would push an atom past its octet; fewer leave a deficit. */
    g_cur_n = 0;
    for (int p = max_pi; p >= 0; p--) {
        enumerate_pi(0, p);
        if (g_cur_n > 0 && g_cur_score[0].excess == 0 && g_cur_score[0].deficit == 0) break;
    }
    if (g_cur_n == 0) return;

    SkeletonScore sk;
    sk.form = g_cur_score[0];
    sk.num_forms = g_cur_n;
    sk.max_heavy_degree = 0;
    sk.h_host_en = 0;
    for (uint8_t k = 0; k < g_n_heavy; k++) {
        uint8_t i = g_heavy[k];
        int degree = g_sigma[i] - g_h_count[k];
        if (degree > sk.max_heavy_degree) sk.max_heavy_degree = degree;
        sk.h_host_en += g_h_count[k] * elem_of(i)->eneg;
    }

    if (g_best_n == 0 || skeleton_better(&sk, &g_best_score)) {
        for (uint8_t k = 0; k < g_cur_n; k++) {
            g_mol->res[k] = g_cur[k];
        }
        g_best_n = g_cur_n;
        g_best_score = sk;
    }
}

static void assign_hydrogens(uint8_t k, uint8_t left)
{
    if (k == g_n_heavy) {
        if (left != 0) return;

        uint8_t next = 0;
        uint8_t pushed = 0;
        bool ok = true;
        for (uint8_t j = 0; j < g_n_heavy && ok; j++) {
            for (uint8_t c = 0; c < g_h_count[j]; c++) {
                if (!push_bond(g_heavy[j], g_hydrogen[next++])) {
                    ok = false;
                    break;
                }
                pushed++;
            }
        }
        if (ok) evaluate_skeleton();
        while (pushed-- > 0) pop_bond();
        return;
    }

    uint8_t room = 0;
    for (uint8_t j = k; j < g_n_heavy; j++) {
        room = (uint8_t)(room + g_h_cap[j]);
    }
    if (room < left) return;

    uint8_t max_c = (g_h_cap[k] < left) ? g_h_cap[k] : left;
    for (int c = max_c; c >= 0; c--) {
        g_h_count[k] = (uint8_t)c;
        assign_hydrogens(k + 1, (uint8_t)(left - c));
    }
    g_h_count[k] = 0;
}

/*
 * Each heavy atom takes at most as many H as its unused natural valence, so
 * acids put H on oxygen. Cations such as NH4+ and H3O+ need one more, so any
 * shortfall is allowed on top of each cap.
 */
static void place_hydrogens(void)
{
    int natural_total = 0;
    for (uint8_t k = 0; k < g_n_heavy; k++) {
        uint8_t i = g_heavy[k];
        int free_valence = natural_valence(g_elem[i]) - g_sigma[i];
        if (free_valence < 0) free_valence = 0;
        g_h_cap[k] = (uint8_t)free_valence;
        natural_total += free_valence;
    }

    int overflow = g_n_hydrogen - natural_total;
    if (overflow < 0) overflow = 0;
    for (uint8_t k = 0; k < g_n_heavy; k++) {
        uint8_t i = g_heavy[k];
        int cap = g_h_cap[k] + overflow;
        int room = max_sigma(g_elem[i]) - g_sigma[i];
        if (cap > room) cap = room;
        g_h_cap[k] = (uint8_t)cap;
        g_h_count[k] = 0;
    }

    assign_hydrogens(0, g_n_hydrogen);
}

/* Attach each non-backbone heavy atom to a backbone atom. */
static void place_terminals(uint8_t k)
{
    if (k == g_n_terminal) {
        place_hydrogens();
        return;
    }

    uint8_t t = g_terminal[k];
    uint8_t start = 0;
    /* Identical terminals are interchangeable: only try non-decreasing hosts. */
    if (k > 0 && g_elem[g_terminal[k - 1]] == g_elem[t]) {
        start = g_term_host[k - 1];
    }
    for (uint8_t h = start; h < g_n_backbone; h++) {
        if (!push_bond(g_backbone[h], t)) continue;
        g_term_host[k] = h;
        place_terminals(k + 1);
        pop_bond();
    }
}

static bool chain_is_mirror_duplicate(void)
{
    for (uint8_t i = 0, j = g_n_backbone - 1; i < j; i++, j--) {
        uint8_t fwd = g_elem[g_backbone[i]];
        uint8_t rev = g_elem[g_backbone[j]];
        if (fwd != rev) return rev < fwd;
    }
    return false;
}

/* Try each distinct element ordering of the backbone chain. */
static void order_backbone(uint8_t pos)
{
    if (pos == g_n_pool) {
        if (chain_is_mirror_duplicate()) return;

        uint8_t made = 0;
        for (uint8_t i = 1; i < g_n_backbone; i++) {
            if (!push_bond(g_backbone[i - 1], g_backbone[i])) break;
            made++;
        }
        if (made + 1 == g_n_backbone) {
            place_terminals(0);
        } else {
            note_failure(INVALID_SKELETON);
        }
        while (made-- > 0) pop_bond();
        return;
    }

    for (uint8_t k = 0; k < g_n_pool; k++) {
        if (g_pool_used[k]) continue;

        bool repeat = false;
        for (uint8_t j = 0; j < k; j++) {
            if (!g_pool_used[j] && g_elem[g_pool[j]] == g_elem[g_pool[k]]) {
                repeat = true;
                break;
            }
        }
        if (repeat) continue;

        g_pool_used[k] = true;
        g_backbone[pos] = g_pool[k];
        order_backbone((uint8_t)(pos + 1));
        g_pool_used[k] = false;
    }
}

/* Stable insertion sort by (natural valence desc, eneg asc) or by element. */
static void sort_atoms(uint8_t *list, uint8_t n, bool by_valence)
{
    for (uint8_t i = 1; i < n; i++) {
        uint8_t v = list[i];
        uint8_t j = i;
        while (j > 0) {
            uint8_t u = list[j - 1];
            bool after;
            if (by_valence) {
                int nu = natural_valence(g_elem[u]);
                int nv = natural_valence(g_elem[v]);
                after = (nu < nv) || (nu == nv && elements[g_elem[u]].eneg > elements[g_elem[v]].eneg);
            } else {
                after = g_elem[u] > g_elem[v];
            }
            if (!after) break;
            list[j] = u;
            j--;
        }
        list[j] = v;
    }
}

static void build_skeletons(void)
{
    g_n_pool = 0;
    g_n_terminal = 0;
    for (uint8_t k = 0; k < g_n_heavy; k++) {
        uint8_t i = g_heavy[k];
        if (is_backbone_elem(g_elem[i])) g_pool[g_n_pool++] = i;
    }

    if (g_n_pool > 0) {
        for (uint8_t k = 0; k < g_n_heavy; k++) {
            uint8_t i = g_heavy[k];
            if (!is_backbone_elem(g_elem[i])) g_terminal[g_n_terminal++] = i;
        }
        sort_atoms(g_pool, g_n_pool, true);
        sort_atoms(g_terminal, g_n_terminal, false);
        memset(g_pool_used, 0, sizeof(g_pool_used));
        g_n_backbone = g_n_pool;
        order_backbone(0);
        return;
    }

    /*
     * No backbone atoms (SO4^2-, ClO4-, O3, Cl2O, XeF4...): one atom is the hub.
     * Candidates are the least electronegative element and any element that
     * appears once; formal charge then decides (Cl2O is O-centered).
     */
    bool has_non_f = false;
    uint8_t min_en = 0xFF;
    for (uint8_t k = 0; k < g_n_heavy; k++) {
        uint8_t elem_idx = g_elem[g_heavy[k]];
        if (elem_idx == ELEM_F) continue;
        has_non_f = true;
        if (elements[elem_idx].eneg < min_en) min_en = elements[elem_idx].eneg;
    }

    for (uint8_t k = 0; k < g_n_heavy; k++) {
        uint8_t hub = g_heavy[k];
        uint8_t elem_idx = g_elem[hub];
        uint8_t count = 0;
        bool first = true;
        for (uint8_t j = 0; j < g_n_heavy; j++) {
            if (g_elem[g_heavy[j]] != elem_idx) continue;
            if (j < k) first = false;
            count++;
        }
        if (!first) continue;
        if (elem_idx == ELEM_F && has_non_f) continue;
        if (has_non_f && elements[elem_idx].eneg != min_en && count != 1) continue;

        g_backbone[0] = hub;
        g_n_backbone = 1;
        g_n_terminal = 0;
        for (uint8_t j = 0; j < g_n_heavy; j++) {
            if (g_heavy[j] != hub) g_terminal[g_n_terminal++] = g_heavy[j];
        }
        sort_atoms(g_terminal, g_n_terminal, false);
        place_terminals(0);
    }
}

/* Display center: most heavy-atom neighbors, then most neighbors, then lowest EN. */
static uint8_t pick_central(const LewisStructure *ls)
{
    int best = -1;
    int best_heavy = -1;
    int best_total = -1;

    for (uint8_t i = 0; i < g_n; i++) {
        if (g_elem[i] == ELEM_H && g_n_heavy > 0) continue;
        int heavy = 0;
        int total = 0;
        for (uint8_t b = 0; b < ls->num_bonds; b++) {
            int other = -1;
            if (ls->bonds[b].a == i) other = ls->bonds[b].b;
            else if (ls->bonds[b].b == i) other = ls->bonds[b].a;
            if (other < 0) continue;
            total++;
            if (g_elem[other] != ELEM_H) heavy++;
        }
        if (best < 0 ||
            heavy > best_heavy ||
            (heavy == best_heavy && total > best_total) ||
            (heavy == best_heavy && total == best_total && elem_of(i)->eneg < elem_of((uint8_t)best)->eneg)) {
            best = i;
            best_heavy = heavy;
            best_total = total;
        }
    }
    return (best < 0) ? 0 : (uint8_t)best;
}

static void generate_single_atom(Molecule *mol)
{
    int shell = is_duet(0) ? 2 : 8;
    if (g_ve != 0 && g_ve != shell) {
        mol->invalid_reason = INVALID_SHELL_RULE;
        return;
    }

    LewisStructure *ls = &mol->res[0];
    memset(ls, 0, sizeof(*ls));
    ls->lone_pairs[0] = (uint8_t)(g_ve / 2);
    ls->formal_charge[0] = (int8_t)(elem_of(0)->valence - g_ve);
    mol->num_res = 1;
}

void generate_resonance(Molecule *mol)
{
    mol->num_res = 0;
    mol->cur_res = 0;
    mol->central = 0;
    mol->invalid_reason = INVALID_NONE;

    if (mol->num_atoms == 0) {
        mol->invalid_reason = INVALID_NO_ATOMS;
        return;
    }

    mol->total_ve = 0;
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        mol->total_ve += elements[mol->atoms[i].elem].valence;
    }
    mol->total_ve -= mol->charge;

    if (mol->total_ve < 0) {
        mol->invalid_reason = INVALID_NEGATIVE_ELECTRONS;
        return;
    }
    if ((mol->total_ve & 1) != 0) {
        mol->invalid_reason = INVALID_ODD_ELECTRONS;
        return;
    }

    g_mol = mol;
    g_n = mol->num_atoms;
    g_ve = mol->total_ve;
    g_n_heavy = 0;
    g_n_hydrogen = 0;
    for (uint8_t i = 0; i < g_n; i++) {
        g_elem[i] = mol->atoms[i].elem;
        if (g_elem[i] == ELEM_H) g_hydrogen[g_n_hydrogen++] = i;
        else g_heavy[g_n_heavy++] = i;
    }

    if (g_n == 1) {
        generate_single_atom(mol);
        return;
    }

    for (uint8_t i = 0; i < g_n; i++) {
        if (max_sigma(g_elem[i]) == 0) {
            mol->invalid_reason = INVALID_SKELETON;
            return;
        }
    }

    /* Lone pairs go to the most electronegative atoms first. */
    for (uint8_t i = 0; i < g_n; i++) {
        uint8_t j = i;
        while (j > 0 && elem_of(g_fill_order[j - 1])->eneg < elem_of(i)->eneg) {
            g_fill_order[j] = g_fill_order[j - 1];
            j--;
        }
        g_fill_order[j] = i;
    }

    memset(&g_ls, 0, sizeof(g_ls));
    memset(g_sigma, 0, sizeof(g_sigma));
    g_best_n = 0;
    g_fail = INVALID_NONE;

    if (g_n_heavy > 0) {
        build_skeletons();
    } else if (g_n == 2 && push_bond(0, 1)) {
        evaluate_skeleton();   /* H2 */
        pop_bond();
    }

    if (g_best_n == 0) {
        mol->invalid_reason = (g_fail == INVALID_NONE) ? INVALID_SKELETON : g_fail;
        return;
    }

    mol->num_res = g_best_n;
    mol->central = pick_central(&mol->res[0]);
}

static void fill_vsepr_fallback(VseprInfo *out)
{
    switch (out->valence_pairs) {
        case 1:
            out->ep_geometry = "Linear";
            out->hybridization = "s";
            out->bond_angle = "180";
            break;
        case 2:
            out->ep_geometry = "Linear";
            out->hybridization = "sp";
            out->bond_angle = "180";
            break;
        case 3:
            out->ep_geometry = "Trigonal Planar";
            out->hybridization = "sp2";
            out->bond_angle = "120";
            break;
        case 4:
            out->ep_geometry = "Tetrahedral";
            out->hybridization = "sp3";
            out->bond_angle = "109.5";
            break;
        case 5:
            out->ep_geometry = "Trigonal Bipyramidal";
            out->hybridization = "sp3d";
            out->bond_angle = "90, 120";
            break;
        case 6:
            out->ep_geometry = "Octahedral";
            out->hybridization = "sp3d2";
            out->bond_angle = "90";
            break;
        case 7:
            out->ep_geometry = "Pentagonal Bipyramidal";
            out->hybridization = "sp3d3";
            out->bond_angle = "72, 90";
            break;
        default:
            out->ep_geometry = "Unknown";
            out->hybridization = "Unknown";
            out->bond_angle = "N/A";
            break;
    }

    if (out->bond_pairs == 0) {
        out->shape = "No Bonded Atoms";
        out->bond_angle = "N/A";
    } else if (out->bond_pairs <= 1) {
        out->shape = "Linear";
        out->bond_angle = "N/A";
    } else {
        out->shape = out->ep_geometry;
    }
}

bool lewis_get_vsepr_info(const Molecule *mol, const LewisStructure *ls, VseprInfo *out)
{
    if (out == NULL) {
        return false;
    }

    memset(out, 0, sizeof(*out));

    if (mol == NULL || ls == NULL) {
        return false;
    }
    if (mol->num_atoms == 0 || mol->central >= mol->num_atoms) {
        return false;
    }

    uint8_t bond_pairs = 0;
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        if (ls->bonds[b].a == mol->central || ls->bonds[b].b == mol->central) {
            bond_pairs++;
        }
    }

    out->bond_pairs = bond_pairs;
    out->lone_pairs = ls->lone_pairs[mol->central];
    out->valence_pairs = (uint8_t)(out->bond_pairs + out->lone_pairs);

    for (uint8_t i = 0; i < (uint8_t)(sizeof(vsepr_rows) / sizeof(vsepr_rows[0])); i++) {
        if (vsepr_rows[i].valence_pairs == out->valence_pairs &&
            vsepr_rows[i].bond_pairs == out->bond_pairs &&
            vsepr_rows[i].lone_pairs == out->lone_pairs) {
            out->ep_geometry = vsepr_rows[i].ep_geometry;
            out->shape = vsepr_rows[i].shape;
            out->hybridization = vsepr_rows[i].hybridization;
            out->bond_angle = vsepr_rows[i].bond_angle;
            return true;
        }
    }

    fill_vsepr_fallback(out);
    return true;
}

const char *invalid_reason_message(InvalidReason reason)
{
    switch (reason) {
        case INVALID_NONE:
            return "No error";
        case INVALID_NO_ATOMS:
            return "No atoms selected";
        case INVALID_NEGATIVE_ELECTRONS:
            return "Invalid charge for selected atoms";
        case INVALID_ODD_ELECTRONS:
            return "Odd electron count (radicals unsupported)";
        case INVALID_SKELETON:
            return "Cannot build a valid bond skeleton";
        case INVALID_LEFTOVER_ELECTRONS:
            return "Could not place all valence electrons";
        case INVALID_SHELL_RULE:
            return "Octet/duet shell constraints failed";
        case INVALID_FORMAL_CHARGE_SUM:
            return "Formal charge sum does not match ion charge";
        default:
            return "No valid structure for this composition";
    }
}
