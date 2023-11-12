## Environment Setup

**Note:** Must use Windows for spectra conversion (Bruker CompassXtract application does not support Mac or Linux)

1. Install (if you haven't already) [Anaconda](https://www.anaconda.com/download/) for **Python** and the **Spyder IDE**
2. Install [ProteoWizard](http://proteowizard.sourceforge.net/)
3. Install **CompassXtract** and **CompassXport** (by [Bruker Daltonics](http://www.bioinfor.com/bruker-data/))



## Running the Conversion Script

1. Open Spyder 
2. Navigate to the directory **"chemcpupy/automation/Bruker"**
3. Open the script **"Bruker_batch_convert.py"**
4. Change the data directory,  `path_source`, to the appropriate address of the data you want to convert 
5. Change the output directory, `path_destin`, to the location you wish to store the converted file
6. Save and run the script



## Example Scripts

To demonstrate how to access converted data using Python or MATLAB, the following scripts have been created:

| Viewing: | Spectra Data                  | Time Data                 |
| -------- | ----------------------------- | ------------------------- |
| Python   | *test_convert_spectra.py*     | *test_convert_raw.py*     |
| MATLAB   | *test_convert_spectra_read.m* | *test_convert_raw_read.m* |



## Conversion Format

The output of the conversions will be an Arch MALDI Mass Spectrometry (**aMALDIms**) file, which is a specifically formatted [HDF5](https://portal.hdfgroup.org/display/HDF5/HDF5) file. HDF5 is well suited for saving and loading large amounts of data, especially when the datasize exceeds a computer's memory cache. 

The structure of the HDF5 output is as follows:

- **\*.hdf5**

  - *attributes* 

    - **filename** : name of the file *(also the name of the original file)*
    - **created** : date file was created *(the conversion file)*
    - **generator** : institution the data was generated at
    - **format** : format name
    - **version** : format version number 
    - **description** : format description
    - **documentation** : link to format documentation (to this document)

  - *groups*
    - **options** : options used for export/conversion
      - *attributes*
        - **include_settings** : exported measurement settings
        - **include_spectra** : exported mass/charge domain data
        - **include_respectra** : exported resampled mass/charge domain data
        - **include_raw** : exported time domain data
        - **compression** : compression used for export (i.e. None, 'gzip', ...)
        - **copy_settings_to_mat** : saved an exhaustive copy of measurement settings to a seperate MAT (MATLAB file)        
        - **save_each_mz_vector** : exported all mass/charge vectors, as opposed to a single one

    - **spectra** : mass/charge-domain data

      - *datasets*
        - **mass_to_charge** : mass/charge array (2D)
        - **signal** : intensity array (2D)
        - **position**: position array
        - **polarity** : ion polarity array (only appears if polarities differ among the spectra)
      - *attributes*
        - **count** : number of measurements
        - **points** : number of points per measurement
        - **created** : date original measurement was taken
        - **ion_polarity** : ion polarity (if equals zero, then polarities differ among the spectra and a polarity dataset has been given instead)  

      - *groups*
        - **parameters**  : spectrum parameters copied from the original '.d' file
        - **statistics** : collection of 1-D arrays (datasets) that are optionally computed upon during conversion to the prevent the need for loading an entire spectrum to get common statistical measures of it
          - *datasets*
            - **max** : maximum instensity of each spectra
            - **min** : minimum instensity of each spectra
            - **mean** : average instensity of each spectra
            - **stdev** : intensity standard deviation of each spectra
            - **median** : 50th intensity percetile of each spectra
            - **quartile1** : 25th intensity percetile of each spectra
            - **quartile3** : 75th intensity percetile of each spectra
            - **sum** : sums of all intensities in a each spectra
            - **area** : area under each spectrum (intensity integrated with trapz over mz)
        - **background** : collections of estimated of background noise levels that are optionally computed during conversion
          - *datasets*
            - **noiseGauss3Sigma** : stand deviation of each spectrum's intensity values that are below the upper limit of the empirical rule
            - **noiseGauss6Sigma** : stand deviation of each spectrum's intensity values that are below twice the upper limit of the empirical rule
            - **noiseTukeyInner** : stand deviation of each spectrum's intensity values that are below the upper limit of the Tukey's inside fences
            - **noiseTukeyOuter** : stand deviation of each spectrum's intensity values that are below the upper limit of the Tukey's outside fences

            - **offsetGauss3Sigma** : average value of each spectrum's intensity values that are below the upper limit of the empirical rule
            - **offsetGauss6Sigma** : average value of each spectrum's intensity values that are below twice the upper limit of the empirical rule
            - **offsetTukeyInner** : average value of each spectrum's intensity values that are below the upper limit of the Tukey's inside fences
            - **offsetTukeyOuter** : average value of each spectrum's intensity values that are below the upper limit of the Tukey's outside fences

            - **cutoffGauss3Sigma** : cutoff *(Mean + 3 Stdev)* below which spectrum points were considered background noise, *noiseGauss3Sigma* is the standard deviation of these points  
            - **cutoffGauss6Sigma** : cutoff *(Mean + 6 Stdev)* below which spectrum points were considered background noise, *noiseGauss6Sigma* is the standard deviation of these points  
            - **cutoffTukeyInner** : cutoff *(Quartile3 + 1.5 InterquartileRange)* below which spectrum points were considered background noise, *noiseTukeyInner* is the standard deviation of these points 
            - **cutoffTukeyOuter** : cutoff *(Quartile3 + 3 InterquartileRange)* below which spectrum points were considered background noise, *noiseTukeyOuter* is the standard deviation of these points
            - **peakcountGauss3Sigma** : number of peaks above the cutoff *cutoffGauss3Sigma*
            - **peakcountGauss6Sigma** : number of peaks above the cutoff *cutoffGauss6Sigma*
            - **peakcountTukeyInner** : number of peaks above the cutoff *cutoffTukeyInner*
            - **peakcountTukeyOuter** : number of peaks above the cutoff *cutoffTukeyOuter*

    - **respectra** : grid resampled mass/charge-domain data

      - *datasets*
        - **mass_to_charge** : mass/charge array (1D)
        - **signal** : resampled intensity array (2D)
        - **position**: position array
        - **polarity** : ion polarity array (only appears if polarities differ among the spectra)
        - **signal_sum** : sum of all intensity vectors (only appears if "collect_sum" is enabled) 

      - *attributes*
        - **count** : number of measurements
        - **points** : number of points per measurement
        - **created** : date original measurement was taken
        - **ion_polarity** : ion polarity (if equals zero, then polarities differ among the spectra and a polarity dataset has been given instead)  
        - **interpolant** : type of interpolation used to syncronize mz vectors
      - *groups*
        - **parameters**  : spectrum parameters copied from the original '.d' file
        - **statistics** : collection of 1-D arrays (datasets) that are optionally computed upon during conversion to the prevent the need for loading an entire spectrum to get common statistical measures of it
          - *datasets*
            - **max** : maximum instensity of each spectra
            - **min** : minimum instensity of each spectra
            - **mean** : average instensity of each spectra
            - **stdev** : intensity standard deviation of each spectra
            - **median** : 50th intensity percetile of each spectra
            - **quartile1** : 25th intensity percetile of each spectra
            - **quartile3** : 75th intensity percetile of each spectra
            - **sum** : sums of all intensities in a each spectra
            - **area** : area under each spectrum (intensity integrated with trapz over mz)
        - **background** : collections of estimated of background noise levels that are optionally computed during conversion
          - *datasets*
            - **noiseGauss3Sigma** : stand deviation of each spectrum's intensity values that are below the upper limit of the empirical rule
            - **noiseGauss6Sigma** : stand deviation of each spectrum's intensity values that are below twice the upper limit of the empirical rule
            - **noiseTukeyInner** : stand deviation of each spectrum's intensity values that are below the upper limit of the Tukey's inside fences
            - **noiseTukeyOuter** : stand deviation of each spectrum's intensity values that are below the upper limit of the Tukey's outside fences

            - **offsetGauss3Sigma** : average value of each spectrum's intensity values that are below the upper limit of the empirical rule
            - **offsetGauss6Sigma** : average value of each spectrum's intensity values that are below twice the upper limit of the empirical rule
            - **offsetTukeyInner** : average value of each spectrum's intensity values that are below the upper limit of the Tukey's inside fences
            - **offsetTukeyOuter** : average value of each spectrum's intensity values that are below the upper limit of the Tukey's outside fences

            - **cutoffGauss3Sigma** : cutoff *(Mean + 3 Stdev)* below which spectrum points were considered background noise, *noiseGauss3Sigma* is the standard deviation of these points  
            - **cutoffGauss6Sigma** : cutoff *(Mean + 6 Stdev)* below which spectrum points were considered background noise, *noiseGauss6Sigma* is the standard deviation of these points  
            - **cutoffTukeyInner** : cutoff *(Quartile3 + 1.5 InterquartileRange)* below which spectrum points were considered background noise, *noiseTukeyInner* is the standard deviation of these points 
            - **cutoffTukeyOuter** : cutoff *(Quartile3 + 3 InterquartileRange)* below which spectrum points were considered background noise, *noiseTukeyOuter* is the standard deviation of these points
            - **peakcountGauss3Sigma** : number of peaks above the cutoff *cutoffGauss3Sigma*
            - **peakcountGauss6Sigma** : number of peaks above the cutoff *cutoffGauss6Sigma*
            - **peakcountTukeyInner** : number of peaks above the cutoff *cutoffTukeyInner*
            - **peakcountTukeyOuter** : number of peaks above the cutoff *cutoffTukeyOuter*
      - *note*: This is an alternative to spectra that has synchonized the mz vectors for all spectra for efficient computation on the 2D signals matrix

    - **raw** : time-domain data
      - *datasets*	
        - **time** : time array
        - **signals** : current array (2D)
        - **position **: position array
      - *attributes*
        - **count** : number of measurements
        - **points** : number of points per measurement
        - **duration** : duration of measurement
        - **sample_rate** : measurement sampling rate
        - **type** : measurement data type
        - **unit_duration** : unit of duration
        - **unit_sample_rate** : unit of rate
    - **scan** : scan information
      - *datasets*
        - **TIC** : total ion current array
        - **elapsedMinutes** : elapsed time array (since start of experiment, in minutes)
        - **position**: position array
        - **number**: scan number array
        - **maxPeak** : max peak intensity array (per spectrum)
    - **method** : method information
      - *groups*
        - **bias** : relevant voltage/bias information
          - *attributes* :
            - **unit** : voltage unit
            - **MALDI_plate_offset** : MALDI plate offset voltage
            - **deflector_plate** : deflector plate voltage
        - **file** : filename information
          - *attributes* :
            - **name** : method name
            - **date** : date method was created
            - **filename** : method filename
            - **filename_original** : original method filename
        - **ion** : ionization information
          - *attributes* :
            - **ionization**: ionization method
            - **polarity** : ion polarity
            - **alternate_polarity** : is ion polarity configured to alternate
        - **laser** : MALDI laser information
          - *attributes* :
            - **power_unit** : laser power unit
            - **power** : laser power
            - **shot_count** : number of laser shots
            - **frequency** : frequency of laser
            - **frequency_unit** : frequency unit
            - **focus** : laser focus
        - **mass** : mass range information
          - *attributes* :
            - **unit** : mass unit
            - **low** : low end of the mass range
            - **high** : high end of the mass range
            - **mode** : mass measurement mode
        - **storage** : acquisition and storage information 
          - *attributes* :
            - **perform_data_reduction** : option to reduce data size
            - **save_full_spectrum** : option to save full spectrum
            - **save_time_series** : option to save time series present
            - **acquisition_size** : acquisition size
            - **processing_size** : processing size
            - **sample_rate** : sample rate of raw data
            - **sample_rate_unit** : sample rate unit
        - **time** : duration information
          - *attributes* :
            - **unit**: time unit
            - **ion_accumulation** : duration ions are accumulated
            - **ion_cooling** : duration ions are cooled
            - **flight_to_detector** : ion time of flight to detector
        - **chroma** : chromatography information
          - *attributes* : 
            - **ESI_high_voltage** : ESI voltage status (on/off)
            - **LC_mode** : LC mode ("MALDI_Automation", "off")
            - **auto_MS_MS** : automatically run MS/MS
            - **source_quench** : quench the source
        - **walk** : MALDI random walk information
          - *attributes* :
            - **enabled** : random walk enabled
            - **pattern** : random walk pattern name
            - **grid_width** : random walk grid width
            - **grid_width_unit** : random walk grid width unit

The created HDF5 file can be viewed in Python using [h5py](https://www.h5py.org/), in MATLAB using [h5read](https://www.mathworks.com/help/matlab/ref/h5read.html), or with the free [HDFVIEW](https://portal.hdfgroup.org/display/HDFVIEW/HDFView) application.



 ## Support

For questions about the converter or format feel free to contact [ChristopherArcadia](https://github.com/ChristopherArcadia).





