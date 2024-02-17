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

For a deeper understanding of the underlying principles and methodologies, you can watch the full lecture here: Efficient Algorithms for Peak Finding - Srini Devadas.
