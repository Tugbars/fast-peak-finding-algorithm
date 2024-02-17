# Peak Finding Algorithm for Impedance Curves
## Overview
This peak finding algorithm is specifically designed for analyzing impedance curves, focusing on identifying the highest peak within a dataset. The algorithm is integral for signal processing and data analysis applications, especially in fields working with impedance measurements where peak characteristics such as prominence and Full Width at Half Maximum (FWHM) are of significant interest.

## Key Features
- High Peak Identification: Determines the highest peak within an impedance curve, crucial for analyzing the material's or circuit's resonant frequency.
- Prominence Calculation: Calculates the prominence of identified peaks to assess their significance against the surrounding data. This helps in distinguishing meaningful peaks from noise or minor fluctuations.
- FWHM Calculation: Computes the Full Width at Half Maximum of the peak, providing insights into the peak's sharpness and the underlying system's damping characteristics.
- Efficient Searching: Utilizes a recursive divide-and-conquer approach to efficiently locate the peak, optimizing performance for large datasets.

## How It Works
- Data Preparation: The algorithm expects data in the form of an array of MqsRawDataPoint_t structures, representing the impedance curve to be analyzed.

- Peak Searching: It employs a recursive method to divide the dataset and pinpoint the highest peak, significantly reducing the search time compared to linear scanning.

- Prominence and FWHM: Once the highest peak is identified, the algorithm calculates its prominence and FWHM to evaluate its physical and analytical relevance.

- Edge Case Handling: Special consideration is given to peaks at the dataset's boundaries, assessing whether a peak might continue beyond the current data range.

## Reference Lecture
The peak finding algorithm implemented in this repository is inspired by the concepts discussed in the lecture by Srini Devadas on efficient algorithms for finding peaks in datasets. The lecture provides a comprehensive overview of the algorithmic approach to identifying peaks within a matrix and its applications in various fields.

For a deeper understanding of the underlying principles and methodologies, you can watch the full lecture here: [Efficient Algorithms for Peak Finding - Srini Devadas](https://youtu.be/HtSuA80QTyo)

## Handling Peak Detection Across Overlapping Arrays
In addition to the primary peak finding algorithm, this repository includes a specialized C method designed to analyze impedance curves across two overlapping arrays, addressing the challenge of capturing peaks that occur at the overlap between these arrays. This method is particularly important for ensuring that no significant peaks are missed due to the segmentation of data across multiple arrays. The algorithm works by considering both arrays simultaneously, thereby enabling the detection of peaks that might not be fully captured within a single array segment.

### Why This Extension Is Necessary
Impedance curves, particularly those with a wide window of data points (e.g., 140 data points), often require a nuanced approach to peak detection. Peaks occurring near the boundary of two arrays may only be partially present in each, necessitating a method that can effectively 'stitch' these segments together for a comprehensive analysis. This extension to the peak finding algorithm precisely addresses this need, ensuring that peaks straddling the boundary between two arrays are accurately identified and analyzed.

### How It Works
The extended algorithm operates by merging the search space across both arrays, effectively treating them as a continuous dataset for the purpose of peak detection. It employs a modified version of the recursive peak finding method to navigate this combined dataset, adjusting its search based on the aggregated data. The algorithm takes into account the following considerations:

Overlap Handling: It dynamically adjusts the search parameters to account for the overlap between arrays, ensuring that the algorithm can seamlessly transition from one array to the next without losing context.

Ignored Indices: Similar to the single-array version, this extended algorithm supports the ignoring of specific indices that have already been evaluated or are deemed irrelevant, thereby optimizing the search process.

Prominence Calculation: For peaks identified at the overlap, the algorithm calculates prominence by considering data points from both arrays, ensuring that the metric accurately reflects the peak's relative prominence within the combined dataset.

FWHM Calculation: The Full Width at Half Maximum (FWHM) for peaks found at the overlap is also calculated across both arrays, providing a true measure of the peak's width irrespective of its position relative to the array boundary.

### Practical Application
This extension is particularly useful in scenarios where data is collected in segments, such as in real-time signal processing or when dealing with large datasets that must be partitioned for analysis. By ensuring that peaks at the boundaries between segments are not overlooked, this approach enhances the accuracy and comprehensiveness of the peak finding process.
