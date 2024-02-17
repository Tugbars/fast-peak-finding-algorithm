#ifndef PEAKPROCESSOR_H
#define PEAKPROCESSOR_H

/*******************************************************************************
 * Includes
 ******************************************************************************/
//#include "mqs/mqs_def.h"

 /*******************************************************************************
  * Defines
  ******************************************************************************/

  /*******************************************************************************
   * Type definitions
   ******************************************************************************/
typedef struct {
	float phaseAngle;
	float impedance;
} MqsRawDataPoint_t;

   /*******************************************************************************
	* Functions
	******************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

	/**
	 * @brief Processes the peak in the given raw data array.
	 *
	 * @param a Pointer to the raw data array.
	 * @param peakIndex Pointer to the variable to store the peak index.
	 * @return true if the peak is successfully processed, false otherwise.
	 */
	bool processPeak(MqsRawDataPoint_t a[], int size, uint16_t *peakIndex, bool* isEdgeCase);

#ifdef __cplusplus
}
#endif

#endif /* PEAKPROCESSOR_H */