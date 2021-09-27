#include "newlisp.h"
#include "protos.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <regex.h>
#include <ctype.h>


char *pte_label[] = {
"X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
"Na", "Mg", "Al", "Si", "P" , "S",  "Cl", "Ar", "K",  "Ca", "Sc",
"Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
"As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
"Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
"Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
"Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
"Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
"Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
"Ds", "Rg"
};

double an2masses[119] = {
0.0,1.00782503207,4.00260325415,7.016004548,9.012182201,11.009305406,
12.0,14.00307400478,15.99491461956,18.998403224,19.99244017542,
22.98976928087,23.985041699,26.981538627,27.97692653246,30.973761629,
31.972070999,34.968852682,39.96238312251,38.963706679,39.962590983,
44.955911909,47.947946281,50.943959507,51.940507472,54.938045141,
55.934937475,58.933195048,57.935342907,62.929597474,63.929142222,
68.925573587,73.921177767,74.921596478,79.916521271,78.918337087,
85.910610729,84.911789737,87.905612124,88.905848295,89.904704416,
92.906378058,97.905408169,98.906254747,101.904349312,102.905504292,
105.903485715,106.90509682,113.90335854,114.903878484,119.902194676,
120.903815686,129.906224399,126.904472681,131.904153457,132.905451932,
137.905247237,138.906353267,139.905438706,140.907652769,141.907723297,
144.912749023,151.919732425,152.921230339,157.924103912,158.925346757,
163.929174751,164.93032207,165.930293061,168.93421325,173.938862089,
174.940771819,179.946549953,180.947995763,183.950931188,186.955753109,
191.96148069,192.96292643,194.964791134,196.966568662,201.970643011,
204.974427541,207.976652071,208.980398734,208.982430435,210.987496271,
222.017577738,222.01755173,228.031070292,227.027752127,232.038055325,
231.03588399,238.050788247,237.048173444,242.058742611,243.06138108,
247.07035354,247.07030708,251.079586788,252.082978512,257.095104724,
258.098431319,255.093241131,260.105504,263.112547,255.107398,259.114500,
262.122892,263.128558,265.136151,281.162061,272.153615,283.171792,283.176451,
285.183698,287.191186,292.199786,291.206564,293.214670};


/*
 * corresponding table of VDW radii.
 * van der Waals radii are taken from A. Bondi,
 * J. Phys. Chem., 68, 441 - 452, 1964,
 * except the value for H, which is taken from R.S. Rowland & R. Taylor,
 * J.Phys.Chem., 100, 7384 - 7391, 1996. Radii that are not available in
 * either of these publications have RvdW = 2.00 Å.
 * The radii for Ions (Na, K, Cl, Ca, Mg, and Cs are based on the CHARMM27
 * Rmin/2 parameters for (SOD, POT, CLA, CAL, MG, CES) by default.
 */
static const double pte_vdw_radius[] = {
    /* X  */ 1.5, 1.2, 1.4, 1.82, 2.0, 2.0,
    /* C  */ 1.7, 1.55, 1.52, 1.47, 1.54,
    /* Na */ 1.36, 1.18, 2.0, 2.1, 1.8,
    /* S  */ 1.8, 2.27, 1.88, 1.76, 1.37, 2.0,
    /* Ti */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Ni */ 1.63, 1.4, 1.39, 1.07, 2.0, 1.85,
    /* Se */ 1.9, 1.85, 2.02, 2.0, 2.0, 2.0,
    /* Zr */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Pd */ 1.63, 1.72, 1.58, 1.93, 2.17, 2.0,
    /* Te */ 2.06, 1.98, 2.16, 2.1, 2.0,
    /* La */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Eu */ 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Er */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* W  */ 2.0, 2.0, 2.0, 2.0, 1.72, 1.66,
    /* Hg */ 1.55, 1.96, 2.02, 2.0, 2.0, 2.0, 2.0,
    /* Fr */ 2.0, 2.0, 2.0, 2.0, 2.0, 1.86,
    /* Np */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Md */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Mt */ 2.0, 2.0, 2.0
};

double get_pte_vdw_radius(int idx)
{
    if ((idx < 1) || (idx >= 112)) return pte_vdw_radius[0];

#if 1
    /* Replace with Hydrogen radius with an "all-atom" radius */
    if (idx == 1)
      return 1.0;    /* H  */
#else
    /* Replace with old VMD atom radii values */
    switch (idx) {
      case  1: return 1.0;    /* H  */
      case  6: return 1.5;    /* C  */
      case  7: return 1.4;    /* N  */
      case  8: return 1.3;    /* O  */
      case  9: return 1.2;    /* F  */
      case 15: return 1.5;    /* P  */
      case 16: return 1.9;    /* S  */
    }
#endif

    return pte_vdw_radius[idx];
}


int get_pte_idx(char *label)
{
  int i;
  char atom[3];
  // zap string
  atom[0] = atom[1] = atom[2] ='\0';

  if(label !=NULL){
    atom[0] = toupper((int) label[0]);
    atom[1] = tolower((int) label[1]);
  }
  if(isdigit(atom[1]) || atom[1]=='_') atom[1] = (char) 0;

  for(i=0; i<112; i++){ // 112 is big enough
    if((pte_label[i][0] == atom[0]) && (pte_label[i][1] == atom[1])) return i;
  }
  return 0;

}

CELL * p_specieNumber (CELL * params){

  char * label;
  size_t labelLen;

  char * labelCopy;
  //CELL * cell;
  int ids = 0;

  getStringSize (params, &label, &labelLen, TRUE);
  if (labelLen>2){
    printf ("Elememt symbol string length > 2 .... may be wrong");
    return (errorProcExt (ERR_STRING_EXPECTED, params->next));
  }
  labelCopy = (char *)allocMemory (labelLen + 1);
  ids = get_pte_idx (label);
  return (stuffInteger (ids));
}

char *get_pte_label(int idx)
{
  if((idx < 1) || (idx >=112)) return pte_label[0];
  return pte_label[idx];
}

CELL * p_num2Specie (CELL * param){
  char *label;
  label = get_pte_label (param->contents);
  return (stuffString (label));
}

double get_pte_mass(int idx)
{
  if((idx < 1) || (idx >=112)) return an2masses[0];
  return an2masses[idx];
}

CELL * p_atommass (CELL *param){
  double mass;
  int idx=0;
  char *label;
  size_t labelLen;
  if(param->type == CELL_STRING){
    getStringSize (param, &label, &labelLen, TRUE);
    idx = get_pte_idx (label);
 }else
    idx = param->contents;

  mass = get_pte_mass (idx);
  return (stuffFloat (mass));
}


CELL * p_atomvdwRadii(CELL * param){
  double radii;
  int idx=0;
  char *label;
  size_t labelLen;
  if(param->type == CELL_STRING){
    getStringSize (param, &label, &labelLen, TRUE);
    idx = get_pte_idx(label);
  }else
      idx = param->contents;
  radii= get_pte_vdw_radius (idx);
  return (stuffFloat(radii));
}


int get_unit_idx(char *label)
{
  int i;
  char atom[3];
  // zap string
  atom[0] = atom[1] = atom[2] ='\0';

  if(label !=NULL){
    atom[0] = toupper((int) label[0]);
    atom[1] = tolower((int) label[1]);
  }
  if(isdigit(atom[1]) || atom[1]=='_') atom[1] = (char) 0;

  for(i=0; i<112; i++){ // 112 is big enough
    if((pte_label[i][0] == atom[0]) && (pte_label[i][1] == atom[1])) return i;
  }
  return 0;
}

double get_unit_value(int idx){
{
    if ((idx < 1) || (idx >= 112)) return pte_vdw_radius[0];

#if 1
    /* Replace with Hydrogen radius with an "all-atom" radius */
    if (idx == 1)
      return 1.0;    /* H  */
#else
    /* Replace with old VMD atom radii values */
    switch (idx) {
      case  1: return 1.0;    /* H  */
      case  6: return 1.5;    /* C  */
      case  7: return 1.4;    /* N  */
      case  8: return 1.3;    /* O  */
      case  9: return 1.2;    /* F  */
      case 15: return 1.5;    /* P  */
      case 16: return 1.9;    /* S  */
    }
#endif

    return pte_vdw_radius[idx];
}
}

CELL * p_unitname(CELL *param){
    double unit;
    char *label;
    int idx = 0;
    size_t labelLen;
    if(param->type == CELL_STRING){
        getStringSize( param, &label, &labelLen, TRUE);
        idx = get_unit_idx(label);
    }
    unit = get_unit_value(idx);
    return (stuffFloat(unit));
}
