#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
//#include <time.h>
//#include <windows.h>
//#include <psapi.h>
#include <stdint.h>

#define MAX_ATTEMPTS 3
#define MAX_IGNORED 3
#define NOISE_TOLERANCE 0.9f
#define PEAK_THRESHOLD 30

static int peakPoint;
static int sweepCounter = 9300;

typedef struct {
	float phaseAngle;
	float impedance;
} MqsRawDataPoint_t;

int recursionCount = 0; // Counter variable

bool shouldBeIgnored(int index, int ignoreIndices[], int numIgnoreIndices)
{
    for (int j = 0; j < numIgnoreIndices; j++)
    {
        if (index == ignoreIndices[j])
        {
            return true;
        }
    }
    return false;
}

float maxrowCombined(MqsRawDataPoint_t a[], int l1, int r1, MqsRawDataPoint_t b[], int l2, int r2, uint16_t *max_index, int *arrayIndex, int ignoreIndices[], int numIgnoreIndices)
{
    float max_val = 0.0f;
    int max_row_index = 0;
    *arrayIndex = 0; // Default to array 'a'

    // Search in array 'a'
    for (int i = l1; i <= r1; i++)
    {
        if (shouldBeIgnored(i, ignoreIndices, numIgnoreIndices))
            continue; // Skip ignored indices
        if (a[i].phaseAngle > max_val)
        {
            max_val = a[i].phaseAngle;
            max_row_index = i;
            *arrayIndex = 1; // Found in array 'a'
        }
    }

    // Search in array 'b'
    for (int i = l2; i <= r2; i++)
    {
        if (shouldBeIgnored(i + r1, ignoreIndices, numIgnoreIndices))
            continue; // Adjust index and skip if needed
        if (b[i].phaseAngle > max_val)
        {
            max_val = b[i].phaseAngle;
            max_row_index = i;
            *arrayIndex = 2; // Found in array 'b'
        }
    }

    *max_index = max_row_index;
    return max_val;
}


static float findPeakrec(MqsRawDataPoint_t a[], int l1, int r1, MqsRawDataPoint_t b[], int l2, int r2, uint16_t *peakIndex, int *arrayIndex, int ignoreIndices[], int numIgnoreIndices)
{
    // Base case for recursion
    if (l1 > r1 && l2 > r2)
    {
        return -1; // No peak found
    }

    float max_val = maxrowCombined(a, l1, r1, b, l2, r2, peakIndex, arrayIndex, ignoreIndices, numIgnoreIndices);

    int mid_combined_a = l1 + (r1 - l1) / 2;
    int mid_combined_b = l2 + (r2 - l2) / 2;

    // Check if the peak is in array 'a'
    if (*arrayIndex == 1 && mid_combined_a > l1 && max_val < a[mid_combined_a - 1].phaseAngle)
    {
        return findPeakrec(a, l1, mid_combined_a - 1, b, l2, r2, peakIndex, arrayIndex, ignoreIndices, numIgnoreIndices);
    }
    // Check if the peak is in array 'b'
    else if (*arrayIndex == 2 && mid_combined_b > l2 && max_val < b[mid_combined_b - 1].phaseAngle)
    {
        return findPeakrec(a, l1, r1, b, l2, mid_combined_b - 1, peakIndex, arrayIndex, ignoreIndices, numIgnoreIndices);
    }
    else
    {
        return max_val; // Peak is found
    }
}

static float calculateProminenceForCombinedArrays(MqsRawDataPoint_t a[], MqsRawDataPoint_t b[], int totalSizeA, int totalSizeB, int arrayIndex, int peakIndex)
{
    float peakValue;

    if (arrayIndex == 1)
    {
        peakValue = a[peakIndex].phaseAngle;
    }
    else if (arrayIndex == 2)
    {
        // Adjust peakIndex for array 'b'
        peakIndex -= totalSizeA;
        peakValue = b[peakIndex].phaseAngle;
    }
    else
    {
        printf("Invalid arrayIndex\n");
        return -1.0; // Return an error value
    }

    float leftMin = peakValue;
    float rightMin = peakValue;

    // Iterate through the left side of the peak
    for (int i = peakIndex - 1; i >= 0; --i)
    {
        float currentValue;
        if (i < totalSizeA)
        { // Use array 'a' if the index is less than totalSizeA
            currentValue = a[i].phaseAngle;
        }
        else
        { // Use array 'b' otherwise
            currentValue = b[i - totalSizeA].phaseAngle;
        }
        if (currentValue < leftMin)
        {
            leftMin = currentValue;
        }
    }
    // Iterate through the right side of the peak
    for (int i = peakIndex + 1; i < totalSizeA + totalSizeB; ++i)
    {
        float currentValue;
        if (i < totalSizeA)
        { // Use array 'a' if the index is less than totalSizeA
            currentValue = a[i].phaseAngle;
        }
        else
        { // Use array 'b' for indices greater than or equal to totalSizeA
            currentValue = b[i - totalSizeA].phaseAngle;
        }
        if (currentValue < rightMin)
        {
            rightMin = currentValue;
        }
    }
    // printf("Inside overlap prominence: leftMin %f and rightMin %f", leftMin, rightMin);

    float prominence = peakValue - fminf(leftMin, rightMin);

    //printf("Prominence: %f\n", prominence);

    return prominence;
}

void calculateSecondOrderDifferenceForCombinedArrays(MqsRawDataPoint_t a[], MqsRawDataPoint_t b[], float secondOrderDiff[], int totalSizeA, int totalSizeB) {
    for (int i = 1; i < totalSizeA + totalSizeB - 1; ++i) {
        float valueA, valueB, valueC;
        if (i < totalSizeA) {
            valueA = a[i].phaseAngle;
            valueC = a[i - 1].phaseAngle;
        } else {
            valueA = b[i - totalSizeA].phaseAngle;
            valueC = b[i - totalSizeA - 1].phaseAngle;
        }
        if (i + 1 < totalSizeA) {
            valueB = a[i + 1].phaseAngle;
        } else {
            valueB = b[i + 1 - totalSizeA].phaseAngle;
        }
        secondOrderDiff[i - 1] = valueB - 2 * valueA + valueC;
    }
}

// Function to find the FWHM peak for combined arrays
static int calculateFWHMForCombinedArrays(MqsRawDataPoint_t a[], MqsRawDataPoint_t b[], int totalSizeA, int totalSizeB, int arrayIndex, int peakIndex, float prominence)
{
    // Calculate the base of the prominence, which is the peak height minus the prominence
    float peakHeight = (arrayIndex == 1) ? a[peakIndex].phaseAngle : b[peakIndex - totalSizeA].phaseAngle;
    float contourLineHeight = peakHeight - prominence;

    // The height at which we measure the FWHM is half the prominence above the contour line
    float halfProminenceHeight = contourLineHeight + (prominence / 2.0f);

    // Initialize left and right indices
    int leftIndex = peakIndex;
    int rightIndex = peakIndex;

    // Iterate through the left side of the peak to find the index where it crosses the half prominence height
    while (leftIndex > 0 && ((arrayIndex == 1 && leftIndex < totalSizeA) ? a[leftIndex].phaseAngle : b[leftIndex - totalSizeA].phaseAngle) > halfProminenceHeight)
    {
        leftIndex--;
    }
    // Adjust index if we are in the second array

    // Iterate through the right side of the peak to find the index where it crosses the half prominence height
    while (rightIndex < totalSizeA + totalSizeB - 1 && ((rightIndex < totalSizeA) ? a[rightIndex].phaseAngle : b[rightIndex - totalSizeA].phaseAngle) > halfProminenceHeight)
    {
        rightIndex++;
    }

    // Calculate FWHM by subtracting indices, considering the contiguous nature of arrays a and b
    int fwhm = fabs(rightIndex - leftIndex);

    return fwhm;
}

static bool isPeakClimbing(MqsRawDataPoint_t b[], int sizeB, int peakIndex, float noiseTolerance)
{
    if (peakIndex <= 0 || peakIndex >= sizeB - 1)
    {
        return false; // Indicating an edge case or error
    }

    int failCount = 0; // Counter for the number of times condition is not met

    for (int i = peakIndex; i < sizeB - 1; i++)
    {
        float derivativeAfter = b[i + 1].phaseAngle - b[i].phaseAngle;

        // Check if the derivative after is less than or equal to the noise tolerance
        if (derivativeAfter <= noiseTolerance)
        {
            failCount++; 
            if (failCount >= 2) 
            {
                return false; // Peak is not climbing if condition failed twice
            }
        }
    }

    // Return true only if failCount is less than 2
    return failCount < 2; 
}

bool processOverlapPeak(MqsRawDataPoint_t *rawData1, int size1, MqsRawDataPoint_t *rawData2, int size2, int maxUpdateAttempts, uint16_t *peakPoint, bool* isEdgeCase)
{
    int peakUpdateAttempts = 0;
    int fwhm = 0;
    uint16_t peakIndex = 0;
    int arrayIndex = -1;
    float peakValue = 0.0f;

    int ignoredIndices[MAX_IGNORED]; // Array to store indices of ignored peaks
    int numIgnored = 0;              // Number of ignored indices

    do
    {
        peakValue = findPeakrec(rawData1, 0, size1 - 1, rawData2, 0, size2 - 1, &peakIndex, &arrayIndex, ignoredIndices, numIgnored);

        peakIndex = (arrayIndex == 1) ? peakIndex : peakIndex + size1;

        // Calculate prominence
        float prominence = calculateProminenceForCombinedArrays(rawData1, rawData2, size1 - 1, size2 - 1, arrayIndex, peakIndex);
        printf("Peak: %f\n", peakValue);
        printf("Index: %d\n", peakIndex);
        printf("p: %f\n", prominence);
        if (prominence > 18.0f)
        {
            int localPeakIndex = arrayIndex == 2 ? peakIndex - size1 : peakIndex;
            // printf("localPeakIndex %d\n", localPeakIndex);
            if (arrayIndex == 2 && (localPeakIndex >= (size2 - PEAK_THRESHOLD)))
            {
                // Adjust peakIndex for array 'b' if necessary
                *isEdgeCase = isPeakClimbing(rawData2, size2, localPeakIndex, NOISE_TOLERANCE);
            }

            fwhm = calculateFWHMForCombinedArrays(rawData1, rawData2, size1, size2, arrayIndex, peakIndex, prominence);
            printf("FWHM: %d\n", fwhm);
            if (fwhm > 15)
            {
                *peakPoint = peakIndex;
                return true;
            }
            else
            {
                printf("FWHM is less than 15.0.\n");

                // Add this peak index to ignored indices
                if (numIgnored < MAX_IGNORED)
                {
                    ignoredIndices[numIgnored++] = peakIndex;
                }

                peakUpdateAttempts++;
                if (peakUpdateAttempts >= maxUpdateAttempts)
                {
                    return false;
                }
            }
        }
        else
        {
            printf("Prominence < 14.0. Not accepting peak.\n");
            return false;
        }

    } while (peakUpdateAttempts < maxUpdateAttempts);
    return false;
}

uint8_t mes_find_overlap_peak(MqsRawDataPoint_t* rawData1, int size1, MqsRawDataPoint_t* rawData2, int size2, int* sweepCounter) {
    uint16_t peakIndex = 0;
    bool isPeakStillClimaxing = false;
    int maxUpdateAttempts = MAX_ATTEMPTS;

    //should return false if isPeakStillClimaxing is true.
    bool peakAccepted = processOverlapPeak(rawData1, size1, rawData2, size2, maxUpdateAttempts, &peakIndex, &isPeakStillClimaxing);
    
    return peakAccepted; // Return the local status determined by the cpocalPeakStatusonditions above
}


int main() {
	float dataset[] = { 10.361000, 10.329520, 10.356401, 10.325025, 10.469888, 10.445896, 10.422787, 10.467480, 10.344401, 10.459909, 10.378614, 10.418076, 10.424760, 10.473890, 10.432741, 10.436613, 10.444571, 10.429080, 10.463049, 10.425678, 10.437474, 10.479097, 10.501722, 10.531240, 10.492681, 10.517651, 10.504417, 10.544653, 10.544653, 10.545215, 10.603968, 10.506781, 10.507369, 10.609545, 10.597960, 10.539934, 10.572769, 10.581369, 10.691141, 10.620659, 10.639743, 10.674317, 10.661292, 10.736961, 10.565084, 10.688236, 10.709663, 10.768684, 10.791526, 10.729278, 10.743296, 10.782402, 10.752879, 10.909691, 10.866303, 10.836424, 10.874863, 10.954317, 10.922943, 10.924746, 10.982296, 10.980767, 10.960667, 11.041705, 10.980650, 10.989566, 11.122129, 11.000278, 11.132257, 11.255452, 11.177774, 11.192039, 11.191874, 11.313030, 11.316112, 11.297583, 11.337660, 11.499168, 11.382261, 11.420565, 11.573527, 11.490598, 11.658082, 11.645509, 11.708488, 11.795426, 11.751255, 11.750044, 11.855704, 11.914387, 12.009725, 11.969546, 12.113441, 12.218554, 12.348103, 12.205872, 12.435554, 12.488775, 12.667537, 12.676172, 12.926952, 12.863553, 12.989057, 13.248148, 13.190160, 13.439136, 13.573619, 13.683957, 13.827342, 13.875702, 14.046788, 14.509664, 14.635375, 14.892009, 14.904869, 15.331629, 15.755693, 15.847921, 16.199364, 16.443979, 16.875294, 17.291578, 17.530399, 18.114887, 18.062302, 18.794970, 19.479204, 19.800901, 21.082626, 20.951014, 22.154087, 22.610720, 23.203785, 24.563568, 25.344297, 26.618078, 27.102108, 28.593575, 29.146513, 30.456078, 31.622009, 32.400932, 34.245522, 35.443687, 36.797287, 37.996586, 38.626411, 39.856213, 40.659065, 41.525280, 41.962757, 42.145386, 41.981716, 41.510342, 41.174747, 40.244114, 38.980572, 37.411938, 36.015099, 34.285168, 32.450775, 30.479216, 28.919357, 28.111219, 27.203331, 25.809673, 25.276243, 23.578642, 22.641386, 21.600714, 21.439640, 20.695690, 19.684826, 19.482126, 18.990290, 17.988312, 18.252808, 17.465487, 16.942823, 16.450624, 16.637707, 16.066063, 15.757387, 15.170953, 15.165143, 14.770429, 14.727147, 14.488015, 14.067205, 13.987227, 13.731712, 13.818885, 13.447730, 13.469353, 13.389613, 13.200713, 13.097751, 12.892175, 13.032427, 12.747318, 12.803812, 12.540964, 12.492415, 12.361678, 12.370881, 12.163138, 12.261773, 11.987444, 11.952088, 11.912817, 11.833737, 12.018749, 11.742359, 11.825325, 11.705390, 11.672668, 11.646121, 11.717649, 11.523814, 11.463550, 11.526981, 11.448123, 11.499317, 11.361500, 11.369127, 11.296580, 11.309932, 11.357458, 11.258648, 11.182965, 11.226593, 11.198554, 11.132441, 11.075950, 11.085775, 11.048738, 11.086349, 11.013202, 11.062451, 10.988196, 10.926581, 10.962508, 10.983298, 11.011072, 10.902027, 10.971194, 10.919538, 10.854755, 10.859086, 10.880175, 10.848403, 10.826693, 10.832817, 10.848177, 10.857426, 10.804535, 10.758336, 10.759258, 10.763223, 10.804464, 10.732544, 10.740483, 10.750152, 10.771185, 10.656355, 10.746325, 10.676956, 10.695798, 10.643116, 10.624805, 10.673359, 10.670972, 10.653358, 10.640178, 10.643605, 10.642442, 10.664634, 10.632175, 10.571341, 10.555463, 10.619086, 10.615108, 10.624764, 10.584524, 10.589610, 10.613992, 10.597569, 10.573765, 10.560243, 10.568216, 10.564842, 10.534982, 10.538974, 10.549685, 10.555965, 10.546945, 10.549246, 10.560552, 10.511511, 10.529139, 10.482478 };

    // Create two arrays with sizes 120 and 180
    MqsRawDataPoint_t rawData1[120];
    MqsRawDataPoint_t rawData2[180];

    // Initialize the arrays with values from the dataset
    for (int i = 0; i < 120; ++i) {
        rawData1[i].phaseAngle = dataset[i];
        rawData1[i].impedance = 0.0;
    }

    for (int i = 0; i < 180; ++i) {
        rawData2[i].phaseAngle = dataset[i + 120];
        rawData2[i].impedance = 0.0;
    }

    //bool peakAccepted = processOverlapPeak(rawData1, 120, rawData2, 180, maxUpdateAttempts, &peakPoint, ignoreIndices, &numIgnoreIndices);
    bool peakAccepted = mes_find_overlap_peak(rawData1, 120, rawData2, 180, &sweepCounter);

    return 0;
}



















