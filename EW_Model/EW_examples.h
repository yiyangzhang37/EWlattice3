#ifndef EW_EXAMPLES_H
#define EW_EXAMPLES_H

#include "EW_Base.h"
#include "EW_BubbleNucl.h"
#include "EW_CSBubble.h"

using namespace Electroweak;


void EW_BaseModel_SymmetricInit(const int n_rows, const int n_cols);

void EW_Nucl_OneBubble(const int n_rows, const int n_cols);

void EW_Nucl_TwoBubbles(const int n_rows, const int n_cols);

void EW_Nucl_NonRand(const int n_rows, const int n_cols);

void EW_Random_Nucl(const int n_rows, const int n_cols);

void EW_Random_Nucl_LongTimeSpectrum(const int n_rows, const int n_cols);

void EW_Random_Nucl_FFT(const int n_rows, const int n_cols);


void EW_CSB_OneBubble(const int n_rows, const int n_cols);


// Efficiency tests
void Lattice_Indexing();

#endif