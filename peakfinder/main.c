/*!
 * Peak Finding Algorithm
 * Author: Tugbars Heptaskin
 * Date: 12/01/2024
 *
 * Description:
 * This set of functions collectively forms a peak finding algorithm designed to identify 
 * and analyze peaks within a dataset. The algorithm is tailored for datasets represented 
 * as arrays of MqsRawDataPoint_t structures.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include "mes_peakfinder.h"


/*!
 * @brief Defines the noise tolerance level for validating edge case climbing peaks.
 *
 * This constant represents the threshold for noise tolerance used in determining whether a peak 
 * is still climbing at the end of a dataset. It is used in the context of peak analysis to 
 * distinguish between genuine rising peaks and minor fluctuations that could be attributed 
 * to noise. A lower value indicates stricter criteria for a peak to be considered as climbing.
 */
#define NOISE_TOLERANCE 0.9f 

/*!
 * @brief Defines the threshold for identifying edge case peaks in a dataset.
 *
 * This constant sets the threshold for defining edge case peaks. An edge case peak is 
 * identified when the peak is near the end of the dataset and has not reached its maximum 
 * (climax) within the dataset's interval. This threshold value determines how close to the 
 * end of the dataset a peak must be to be considered an edge case. It is used to decide 
 * whether to check if a peak is still climbing or if it may continue in a subsequent dataset.
 */
#define PEAK_THRESHOLD  30 



/*!
 * @brief Calculates the prominence of a peak in a dataset.
 *
 * This function computes the prominence of a specified peak in an array of data points. 
 * Prominence in this context refers to the height of the peak relative to the lowest contour 
 * line that encloses the peak and no higher peak. It is a measure of how a peak stands out 
 * from the surrounding baseline and is important in signal processing and data analysis 
 * for distinguishing significant peaks from minor fluctuations.
 *
 * The function first identifies the nearest higher peaks (or the ends of the dataset if no 
 * higher peaks are present) on both the left and right sides of the specified peak. It then 
 * finds the minimum value within this range, which represents the base of the peak. The 
 * prominence is calculated as the difference between the peak's value and this minimum value.
 *
 * @param a The array of data points (MqsRawDataPoint_t) in which the peak is located.
 * @param size The size of the array.
 * @param peakIndex The index of the peak within the array.
 * @return The prominence of the specified peak.
 */
static float findProminence(MqsRawDataPoint_t a[], int size, int peakIndex)
{
    // Initialize variables to track the nearest higher peaks or ends
    int leftBoundary = 0;
    int rightBoundary = size - 1;

    float peak_val = a[peakIndex].phaseAngle;

    // Find the nearest higher peak or end on the left
    for (int i = peakIndex - 1; i >= 0; i--)
    {
        if (a[i].phaseAngle > peak_val)
        {
            leftBoundary = i;
            break;
        }
    }

    // Find the nearest higher peak or end on the right
    for (int i = peakIndex + 1; i < size; i++)
    {
        if (a[i].phaseAngle > peak_val)
        {
            rightBoundary = i;
            break;
        }
    }

    // Find the minimum value within the boundaries
    float minValue = a[rightBoundary].phaseAngle;
    for (int i = leftBoundary; i <= rightBoundary; i++)
    {
        if (a[i].phaseAngle < minValue)
        {
            minValue = a[i].phaseAngle;
        }
    }

    // printf("min Value %f", minValue);

    // Calculate and return the prominence
    return peak_val - minValue;
}

/*!
 * @brief Finds the index of the maximum value in a column of a 2D array, ignoring specified indices.
 *
 * @param a The array of data points (MqsRawDataPoint_t) to search through.
 * @param size The number of elements in the array.
 * @param col The column in the array to search for the maximum value.
 * @param max_val A pointer to store the maximum value found.
 * @param max_index A pointer to store the index of the maximum value.
 * @param ignoreIndices An array of indices to be ignored during the search.
 * @param numIgnoreIndices The number of indices to ignore.
 * @return The index of the maximum value found in the specified column.
 */
static int maxrow(MqsRawDataPoint_t a[], int size, int col, float *max_val, int *max_index, int ignoreIndices[], int numIgnoreIndices)
{
    for (int i = 0; i < size; i++)
    {
        // Skip the ignored indices
        int ignore = 0;
        for (int j = 0; j < numIgnoreIndices; j++)
        {
            if (i == ignoreIndices[j])
            {
                ignore = 1;
                break;
            }
        }

        if (ignore)
        {
            continue;
        }

        if (*max_val < a[i].phaseAngle)
        {
            *max_val = a[i].phaseAngle;
            *max_index = i;
        }
    }
    return *max_index;
}

/*!
 * @brief Recursively finds a peak in a dataset using a divide-and-conquer approach.
 *
 * This function implements a recursive peak finding algorithm. It divides the dataset 
 * into two halves at each recursive step and determines the direction (left or right) 
 * to continue the search based on the comparison of adjacent values. This divide-and-conquer 
 * approach significantly reduces the time complexity compared to a linear search, improving 
 * performance, especially in large datasets.
 *
 * The function also supports ignoring specific indices in the dataset, which can be useful 
 * in cases where certain data points have low FWHM. 
 *
 * @param a The array of data points (MqsRawDataPoint_t) to search through for a peak.
 * @param size The size of the array.
 * @param l The starting index of the current search window.
 * @param r The ending index of the current search window.
 * @param peakIndex A pointer to store the index of the found peak.
 * @param ignoreIndices An array of indices to be ignored during the search.
 * @param numIgnoreIndices The number of indices to ignore.
 * @return The value of the peak found, or -1 if no peak is found.
 */
static float findPeakRec(MqsRawDataPoint_t a[], int size, int l, int r, uint16_t *peakIndex, int ignoreIndices[], int numIgnoreIndices)
{

    if (l > r)
        return -1;

    int mid = (l + r) / 2;
    float max_val = 0.0f;
    int max_index = 0;

    // Skip the ignored indices in the maxrow function
    int max_row_index = maxrow(a, size, mid, &max_val, &max_index, ignoreIndices, numIgnoreIndices);

    // printf("%f ", a[max_row_index].phaseAngle);

    if (mid == 0 || mid == size - 1)
    {
        *peakIndex = max_row_index;
        return max_val;
    }

    if (max_val < a[mid - 1].phaseAngle)
        return findPeakRec(a, size, l, mid - 1, peakIndex, ignoreIndices, numIgnoreIndices);
    else if (max_val < a[mid + 1].phaseAngle)
        return findPeakRec(a, size, mid + 1, r, peakIndex, ignoreIndices, numIgnoreIndices);
    else
    {
        *peakIndex = max_row_index;
        return max_val;
    }
}

/*!
 * @brief Calculates the Full Width at Half Maximum (FWHM) of a peak in a dataset.
 *
 * This function computes the Full Width at Half Maximum (FWHM) for a specified peak. 
 * FWHM is a measure of the width of a peak at its half maximum height. This implementation 
 * calculates the FWHM based on half the prominence of the peak, using a methodology similar 
 * to that described in MathWorks' findpeaks function documentation.
 *
 * The process involves finding the height at half the prominence above the contour line
 * (base level) of the peak. The function then locates the left and right indices where 
 * the signal crosses this half-prominence height. The FWHM is determined as the distance 
 * between these two indices.
 *
 * Special consideration is given in cases where the signal does not exactly cross the 
 * half-prominence height at a data point, by interpolating between points if necessary.
 * for more information: // https://www.mathworks.com/help/signal/ref/findpeaks.html#buhd6xj
 *
 * @param a The array of data points (MqsRawDataPoint_t) containing the peak.
 * @param size The size of the array.
 * @param peakIndex The index of the peak within the array.
 * @param prominence The prominence of the peak, used to determine the half-prominence height.
 * @return The FWHM of the specified peak, calculated based on half the prominence.
 */
static int calculateFWHM(MqsRawDataPoint_t a[], int size, int peakIndex, float prominence)
{
    // First, find the base of the peak
    float peakHeight = a[peakIndex].phaseAngle;
    float contourLineHeight = peakHeight - prominence;

    // The height at which we measure the FWHM is half the prominence above the contour line
    float halfProminenceHeight = contourLineHeight + (prominence / 2.0f);

    // Find the left and right indices where the phase angle crosses the half-prominence height
    int leftIndex = peakIndex;
    while (leftIndex > 0 && a[leftIndex].phaseAngle > halfProminenceHeight)
    {
        leftIndex--;
    }

    int rightIndex = peakIndex;
    while (rightIndex < size - 1 && a[rightIndex].phaseAngle > halfProminenceHeight)
    {
        rightIndex++;
    }

    // Calculate FWHM using the phase angles at left and right indices
    // Be sure to handle the cases where the signal doesn't cross the half-prominence height
    // at a data point exactly by interpolating between points if necessary
    int fwhm = fabsf(rightIndex - leftIndex);

    return fwhm;
}

/*!
 * @brief Determines if a peak is still climbing at the end of a dataset.
 *
 * This function assesses whether the identified peak in a dataset is still rising 
 * as it reaches the end of the dataset. This is important in peak finding algorithms,
 * particularly when analyzing segments of data where the peak might extend beyond 
 * the current dataset's boundary.
 *
 * The function iterates from the peak index to the end of the dataset, calculating 
 * the derivative (rate of change) at each point. It checks if this derivative is 
 * less than or equal to a specified noise tolerance. If the condition fails more than once,
 * it indicates that the peak is no longer climbing.
 *
 * This check helps to determine if the current dataset's peak is part of a larger peak 
 * that might be fully realized in subsequent datasets. If the peak is still climbing 
 * at the end of the current dataset, there may be a need to analyze the next dataset 
 * to find the true peak.
 *
 * @param b The array of data points (MqsRawDataPoint_t) containing the peak.
 * @param sizeB The size of the array.
 * @param peakIndex The index of the peak within the array.
 * @param noiseTolerance The tolerance level for the derivative to be considered noise.
 * @return True if the peak is still climbing; false otherwise.
 */
static bool isPeakClimbing(MqsRawDataPoint_t b[], int sizeB, int peakIndex, float noiseTolerance)
{
    if (peakIndex <= 0 || peakIndex >= sizeB - 1)
    {
        return false; 
    }

    int failCount = 0; // Counter for the number of times condition is not met

    for (int i = peakIndex; i < sizeB - 1; i++)
    {
        float derivativeAfter = b[i + 1].phaseAngle - b[i].phaseAngle;

        // Check if the derivative after is less than or equal to the noise tolerance
        if (derivativeAfter <= noiseTolerance)
        {
            failCount++; 
            if (failCount >= 2) // Check if it's the second time
            {
                return false; // Peak is not climbing if condition failed twice
            }
        }
    }

    // Return true only if failCount is less than 2
    return failCount < 2; 
}

float calculateDampingRatio(float resonanceFrequency, float FWHM) {
    float dampingRatio = resonanceFrequency / (2 * M_PI * FWHM);
    return dampingRatio;
}

double lorentzian(double frequency, double peakHeight, double resonanceFrequency, double halfWidthAtHalfMaximum) {
    return (peakHeight / M_PI) * (halfWidthAtHalfMaximum / (pow(frequency - resonanceFrequency, 2) + pow(halfWidthAtHalfMaximum, 2)));
}

/*!
 * @brief Processes and validates a peak within a dataset.
 *
 * This function identifies and validates a peak in a given dataset. The peak is first identified
 * using a recursive peak-finding algorithm. Once found, the function calculates the peak's
 * prominence and Full Width at Half Maximum (FWHM) to determine its significance and breadth.
 *
 * The peak is considered valid if:
 *   - The prominence exceeds a specified threshold, indicating it is a significant peak.
 *   - The FWHM is greater than a certain value, ensuring the peak is not too narrow.
 *
 * Additionally, if the peak is near the end of the dataset, the function checks if the peak is 
 * still climbing, indicating that it might continue in the next dataset. This is determined using 
 * the `isPeakClimbing` function.
 *
 * If the peak does not meet these criteria, it is skipped, and the function attempts to find 
 * another peak, up to a maximum number of attempts. Peaks that are skipped are recorded in an 
 * array to prevent reprocessing in subsequent attempts.
 *
 * @param a The array of data points (MqsRawDataPoint_t) containing the potential peak.
 * @param size The size of the array.
 * @param peakIndex A pointer to store the index of the identified peak.
 * @param isEdgeCase A pointer to a boolean flag indicating if the peak is an edge case.
 * @return True if a valid peak is found and processed; false otherwise.
 */
bool processPeak(MqsRawDataPoint_t a[], int size, uint16_t *peakIndex, bool* isEdgeCase)
{
    int skippedIndices[3]; // Array to store the indices of skipped peaks
    int skippedCount = 0;  // Count of skipped peaks
    int maxAttempts = 3;   // Maximum number of attempts
    int fwhm = 0;
    int retry = 0;

    do
    {
        float peakValue = findPeakRec(a, size, 0, size - 1, peakIndex, skippedIndices, skippedCount);

        if (peakValue == -1)
        {
            printf("No peak found.\n");
            return false;
        }

        printf("\nPeak: %f\n", peakValue);
        printf("Index: %d\n", *peakIndex);

        // Check prominence
        float prominence = findProminence(a, size - 1, *peakIndex);
        printf("Prominence: %f\n", prominence);

        if (prominence > 18.0f)
        {
            // Check FWHM
            fwhm = calculateFWHM(a, size, *peakIndex, prominence);
            printf("FWHM: %d\n", fwhm);

            // Check if peak is near the end and potentially still climaxing
            if (*peakIndex >= size - PEAK_THRESHOLD)
            {
                *isEdgeCase = isPeakClimbing(a, size, *peakIndex, NOISE_TOLERANCE);
            }

            if (fwhm > 15) 
            {
                return true; // Peak accepted
            }
            else
            {
                printf("FWHM is less than 15.0. Retrying...\n");
                // Store the index of the skipped peak
                if (skippedCount < 3)
                {
                    skippedIndices[skippedCount++] = *peakIndex;
                }
            }
        }
        else
        {
            printf("Prominence < 18.0.\n");
            // Exit the loop if the prominence is less than 14.0
            break;
        }

        retry++;
    } while (retry < maxAttempts);

    return false;
}

bool mes_find_peak(MqsRawDataPoint_t* rawData, int size, int* sweepCounter) {
    uint16_t peakIndex = 0;
    bool isPeakStillClimaxing = false;
   
    bool peakAccepted = processPeak(rawData, size, &peakIndex, &isPeakStillClimaxing);
    
    return peakAccepted; // Return the updated status value.
}

int main() {
    float dataset[301] = { 10.361000, 10.329520, 10.356401, 10.325025, 10.469888, 10.445896, 10.422787, 10.467480, 10.344401, 10.459909, 10.378614, 10.418076, 10.424760, 10.473890, 10.432741, 10.436613, 10.444571, 10.429080, 10.463049, 10.425678, 10.437474, 10.479097, 10.501722, 10.531240, 10.492681, 10.517651, 10.504417, 10.544653, 10.544653, 10.545215, 10.603968, 10.506781, 10.507369, 10.609545, 10.597960, 10.539934, 10.572769, 10.581369, 10.691141, 10.620659, 10.639743, 10.674317, 10.661292, 10.736961, 10.565084, 10.688236, 10.709663, 10.768684, 10.791526, 10.729278, 10.743296, 10.782402, 10.752879, 10.909691, 10.866303, 10.836424, 10.874863, 10.954317, 10.922943, 10.924746, 10.982296, 10.980767, 10.960667, 11.041705, 10.980650, 10.989566, 11.122129, 11.000278, 11.132257, 11.255452, 11.177774, 11.192039, 11.191874, 11.313030, 11.316112, 11.297583, 11.337660, 11.499168, 11.382261, 11.420565, 11.573527, 11.490598, 11.658082, 11.645509, 11.708488, 11.795426, 11.751255, 11.750044, 11.855704, 11.914387, 12.009725, 11.969546, 12.113441, 12.218554, 12.348103, 12.205872, 12.435554, 12.488775, 12.667537, 12.676172, 12.926952, 12.863553, 12.989057, 13.248148, 13.190160, 13.439136, 13.573619, 13.683957, 13.827342, 13.875702, 14.046788, 14.509664, 14.635375, 14.892009, 14.904869, 15.331629, 15.755693, 15.847921, 16.199364, 16.443979, 16.875294, 17.291578, 17.530399, 18.114887, 18.062302, 18.794970, 19.479204, 19.800901, 21.082626, 20.951014, 22.154087, 22.610720, 23.203785, 24.563568, 25.344297, 26.618078, 27.102108, 28.593575, 29.146513, 30.456078, 31.622009, 32.400932, 34.245522, 35.443687, 36.797287, 37.996586, 38.626411, 39.856213, 40.659065, 41.525280, 41.962757, 42.145386, 41.981716, 41.510342, 41.174747, 40.244114, 38.980572, 37.411938, 36.015099, 34.285168, 32.450775, 30.479216, 28.919357, 28.111219, 27.203331, 25.809673, 25.276243, 23.578642, 22.641386, 21.600714, 21.439640, 20.695690, 19.684826, 19.482126, 18.990290, 17.988312, 18.252808, 17.465487, 16.942823, 16.450624, 16.637707, 16.066063, 15.757387, 15.170953, 15.165143, 14.770429, 14.727147, 14.488015, 14.067205, 13.987227, 13.731712, 13.818885, 13.447730, 13.469353, 13.389613, 13.200713, 13.097751, 12.892175, 13.032427, 12.747318, 12.803812, 12.540964, 12.492415, 12.361678, 12.370881, 12.163138, 12.261773, 11.987444, 11.952088, 11.912817, 11.833737, 12.018749, 11.742359, 11.825325, 11.705390, 11.672668, 11.646121, 11.717649, 11.523814, 11.463550, 11.526981, 11.448123, 11.499317, 11.361500, 11.369127, 11.296580, 11.309932, 11.357458, 11.258648, 11.182965, 11.226593, 11.198554, 11.132441, 11.075950, 11.085775, 11.048738, 11.086349, 11.013202, 11.062451, 10.988196, 10.926581, 10.962508, 10.983298, 11.011072, 10.902027, 10.971194, 10.919538, 10.854755, 10.859086, 10.880175, 10.848403, 10.826693, 10.832817, 10.848177, 10.857426, 10.804535, 10.758336, 10.759258, 10.763223, 10.804464, 10.732544, 10.740483, 10.750152, 10.771185, 10.656355, 10.746325, 10.676956, 10.695798, 10.643116, 10.624805, 10.673359, 10.670972, 10.653358, 10.640178, 10.643605, 10.642442, 10.664634, 10.632175, 10.571341, 10.555463, 10.619086, 10.615108, 10.624764, 10.584524, 10.589610, 10.613992, 10.597569, 10.573765, 10.560243, 10.568216, 10.564842, 10.534982, 10.538974, 10.549685, 10.555965, 10.546945, 10.549246, 10.560552, 10.511511, 10.529139, 10.482478 };
    int sweepCounter = 9300;
    
    // Create an array of MqsRawDataPoint_t with phaseAngle values from the dataset
    MqsRawDataPoint_t rawData[301];
    for (int i = 0; i < 301; ++i) {
        rawData[i].phaseAngle = dataset[i];
        rawData[i].impedance = 0.0;  // You can set the impedance to a default value
    }

    //bool peakAccepted = processPeak(rawData, 301, &peakPoint, ignoreIndices, &numIgnoreIndices);
    bool peakAccepted = mes_find_peak(rawData, 301, &sweepCounter);
 
    return 0;
}















