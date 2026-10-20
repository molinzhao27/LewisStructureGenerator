#ifndef LEWIS_MODEL_H
#define LEWIS_MODEL_H

#include <stdbool.h>
#include <stdint.h>

/* Screen layout constants */
#define SCR_W           320
#define SCR_H           240

/* Top info bar */
#define INFO_Y          0
#define INFO_H          30

/* Selected atoms bar */
#define SEL_Y           32
#define SEL_H           43

/* Periodic table area */
#define PT_Y            80
#define PT_CELL_W       17
#define PT_CELL_H       16
#define PT_COLS         18
#define PT_ROWS         5
#define PT_LEFT_PAD     2

/* Lewis structure screen */
#define LEWIS_CENTER_X  160
#define LEWIS_CENTER_Y  135
#define BOND_LEN        55
#define DOT_R           2
#define DOT_DIST        14

/* Molecule limits */
#define MAX_ATOMS       12
#define MAX_HEAVY       6
#define MAX_BONDS       12
#define MAX_RESONANCE   6

/* Frame rate target */
#define TARGET_FPS      20
#define FRAME_TICKS     (32768 / TARGET_FPS)

/* Color palette indices (default TI palette) */
#define COL_BLACK       0x00
#define COL_WHITE       0xFF
#define COL_RED         0xE0
#define COL_BLUE        0x10
#define COL_GREEN       0x04
#define COL_YELLOW      0xE7
#define COL_CYAN        0x1F
#define COL_MAGENTA     0xE3
#define COL_ORANGE      0xE4
#define COL_GRAY        0xB5
#define COL_DKGRAY      0x6B

/* Element info card constants */
#define CARD_X          56
#define CARD_Y          84
#define CARD_W          86
#define CARD_H          80

/* Index constants for quick lookup */
#define ELEM_NONE       0xFF
#define NUM_ELEMENTS    34

/* Stable element index aliases (match elements[] order) */
#define ELEM_H          0
#define ELEM_HE         1
#define ELEM_C          5
#define ELEM_N          6
#define ELEM_O          7
#define ELEM_F          8
#define ELEM_S          15

typedef struct {
    char     symbol[3];
    char     name[14];
    uint8_t  atomic_num;   /* atomic number (Z) */
    uint8_t  valence;      /* valence electrons */
    uint8_t  bond_cap;     /* typical bonding capacity */
    uint8_t  eneg;         /* electronegativity * 10 (integer) */
    uint8_t  period;       /* 1-based */
    uint8_t  group;        /* 1-based, 1-18 */
    uint8_t  color;        /* display color */
} Element;

extern const Element elements[NUM_ELEMENTS];
extern uint8_t pt_grid[PT_ROWS][PT_COLS];

void init_pt_grid(void);

typedef struct {
    uint8_t  elem;         /* index into elements[] */
} Atom;

typedef struct {
    uint8_t  a, b;         /* atom indices */
    uint8_t  order;        /* 1=single, 2=double, 3=triple */
} Bond;

typedef struct {
    uint8_t  lone_pairs[MAX_ATOMS]; /* lone pairs per atom */
    Bond     bonds[MAX_BONDS];
    uint8_t  num_bonds;
    int8_t   formal_charge[MAX_ATOMS];
} LewisStructure;

typedef enum {
    INVALID_NONE = 0,
    INVALID_NO_ATOMS,
    INVALID_NEGATIVE_ELECTRONS,
    INVALID_ODD_ELECTRONS,
    INVALID_SKELETON,
    INVALID_LEFTOVER_ELECTRONS,
    INVALID_SHELL_RULE,
    INVALID_FORMAL_CHARGE_SUM
} InvalidReason;

typedef struct {
    Atom     atoms[MAX_ATOMS];
    uint8_t  num_atoms;
    int8_t   charge;       /* overall molecular charge */

    /* Generated structures */
    LewisStructure res[MAX_RESONANCE];
    uint8_t  num_res;
    uint8_t  cur_res;      /* currently displayed resonance form */

    uint8_t  central;      /* index of central atom */
    int      total_ve;     /* total valence electrons */
    InvalidReason invalid_reason;
} Molecule;

void molecule_reset(Molecule *mol);

#endif
