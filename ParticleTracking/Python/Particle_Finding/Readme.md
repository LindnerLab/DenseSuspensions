# Particle Finding Algorithm
This is the documentation for a particle finding algorithm. It is based on a
the convolution of the image with particles to be detected and a template of
a particle to be detected (mask). For now, this algorithm only works with
rotationally symmetric particles (circles/annuli).

All functions in *particle_finding_functions* are stand alone. However, as
input is expected in a certain format (explained in detail in
*particle_finding_script*, and in the docstrings of the functions themselves),
it is recommended to run it from the provided example script and adapt the
parameters where needed.

Furthermore, if one desired the run the particle finding algorithm in parallel
and required image pretreatment functions not present implemented. It is
recommended to copy the *find_parallel* function (including the @ray.remote
environment) and change the *img_pretreatment* function to the desired
pretreatment (the output should be a normalized image).

## Walkthough of example script
- First, the desired folder structure has to be provided
- Next, the script builds a dict with the desired settings, or it loads a
  previously saved settings dict.
- Now, one of three processes is started:
    1. If **parallel_processing == True**, all files will be passed to
    *find_parallel*, and will be processes in parallel on *settings\['nCores'\]*
    cores.
    1. If **verbose == True**, the files specified by *files_of_interest* will
    be passed to *image_pretreatment* and then to *find_serial* and will be
    processed one by one. Besides saving the particle centroids found, the
    centroids will also be returned in *particles*, which is a list of a list
    of numpy arrays, where the index in the top list indicates the image used,
    the index in the sublist indicates the particle type found (given
    *settings\['R'\]*), and the first column of the numpy array contains the
    x-coordinates of the particles, and the second column the y-coordinates.
    All coordinates are in pixels, with the origin at the top-left corner of
    the image.
    1. Or if both are False, all images are passed to *find_serial*. Output
    is only saved if specified in settings*\['save_files'\]* (not returned as
    list)
- The particle finding works as follows:
    1. First, the normalized image is thresholded (where all values below the
    threshold are set to 0, and values above are kept as is), using
    *settings\['thresh_img'\]*.
    1. Then the appropriate mask is made. The size of the entire mask is
    identical to the image size, and the size of the ROI is indicated by
    *settings\['R'\]*.
    1. After that, the image and mask are FFT'ed, then multiplied, followed
    by an IFFT of the convoluted image.
    1. Following the convolution, the convoluted image is binarized, using
    *settings\['thresh_conv'\]* as a threshold
    1. Next, possible particle locagittions are labelled using scipy.measure
    regeionprops.
    1. ROI's are then separated from false positives given
    *settings\['selection_criteria'\]*.
    1. The centroid of the ROI's is then determined using regionprops
    weighted_centroid, where each pixel in the ROI is weighted according to
    intensity of the convoluted image.
    1. These centroids are then saved/returned as list/discarded, depending on
    the settings given.
    
# Contributors
**Contributors** : Lars Kool

**Affiliations** : Laboratoire Physique et Méchanique des Milieux Hétérogènes
(PMMH), ESPCI, Paris, France

# Funding Acknowledgement
This project has received funding from the European Union’s Horizon 2020
research and innovation programme under the Marie Skłodowska-Curie grant
agreement No. 813162

More info on this programme: https://caliper-itn.org/
<img src="https://caliperitn.files.wordpress.com/2019/06/cropped-a-13-2.png height=200">
<img src="https://eacea.ec.europa.eu/sites/eacea-site/files/flag_2colors.png" height=200>
