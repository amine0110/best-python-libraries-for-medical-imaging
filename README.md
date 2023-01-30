# The best Python Libraries for Medical Imaging
## Pydicom, Nibabel, VTK,…


In this blog post, we will discuss the best libraries that can be used in Python for medical imaging.
Original blog post on my website [here](https://pycad.co/the-best-python-libraries-for-medical-imaging/).

Learn how to effectively manage and process DICOM files in Python with our comprehensive course, designed to equip you with the skills and knowledge you need to succeed.

https://www.learn.pycad.co/course/dicom-simplified

# Outline

- [Pydicom](https://github.com/amine0110/best-python-libraries-for-medical-imaging#pydicom)
- [Nibabel](https://github.com/amine0110/best-python-libraries-for-medical-imaging#nibabel)
- [Dicom2nifti](https://github.com/amine0110/best-python-libraries-for-medical-imaging#dicom2nifti)
- [SimpleITK](https://github.com/amine0110/best-python-libraries-for-medical-imaging#simpleitk)
- [VTK](https://github.com/amine0110/best-python-libraries-for-medical-imaging#vtk)
- [Numpy to STL (using NumPy-stl)](https://github.com/amine0110/best-python-libraries-for-medical-imaging#numpy-stl)
- [MedPy](https://github.com/amine0110/best-python-libraries-for-medical-imaging#medpy)
- [MONAI](https://github.com/amine0110/best-python-libraries-for-medical-imaging#monai)

---

# Pydicom

[🔗](https://pydicom.github.io/pydicom/stable/old/base_element.html#tag) Pydicom is an open-source library for working with Dicom files. It's what we use to load, edit, and save Dicom files!

Examples:

```Python
import pydicom
# 1. Load a dicom file
ds = pydicom.dcmread('path to dicom file')

# 2. Extract the patient's name
ds.PatientName

# 3. Extract the patient's birthday
ds.PatientBirthDate

# 4. Extract the study ID
ds.StudyID

# 5. Extract the study date
ds.StudyDate

# 6. Extract the image from the dicom (image to display)
ds.pixel_array
```


If we want to show/visualize the array extracted from the Dicom file, we can use the following code:

```Python
import pydicom
from PIL import Image
import numpy as np

ds = pydicom.dcmread('the path to the dicom file')

# Convert to fload to avoid overflow or underflow losses
image_2d = ds.pixel_array.astype(float)

# Rescaling grey scale between 0-255
image_2d_scaled = (np.maximum(image_2d, 0) / (image_2d.max()) * 255.0

# convert to Uint
umage_2d_scaled = np.uint8(image_2d_scaled)

im = Image.fromarray(image_2d_scaled)

# Show the image
im.show()

# Save the image as PNG or JPG
im.save('patient_slice.png') # or .jpg
```

---

# Nibabel

[🔗](https://nipy.org/nibabel/) Then there's nibabel, which is the most commonly used library for dealing with nifti files. There are many options in this library; here are a few to get you started:

```Python
import nibale as nib

# 1. Load/upload a nifti file
nifti_file = nib.load('The path to the nifti file')

# 2. Return the header of the file (which contains informations about the voxels, image dimensions...)
header = nifti_file.header

# 3. Return the 3D array containing slices
array_3d = nifti_file.get_fdata()

# 4. Get the affine matrix
nifti_file.affine

# 5. Create a nifti file from numpy array
new_nifti_file = nib.nifti1.Nifti1Image(data, affine, header=new_header)

# 6. Save the new nifti file
save(new_nifti_fle, 'new_nifti_file.nii') # or .nii.gz

```

---

# Dicom2nifti

[🔗](https://pypi.org/project/dicom2nifti/) This is a very useful library for converting Dicom series into nifti files with a single function.

```Python
import dicom2nifti

dicom2nifti.dicom_series_to_nifti('The path to the dicom series', 'output_nifti_file.nii') # or .nii.gz

```

---

# SimpleITK

[🔗](https://simpleitk.readthedocs.io/en/master/) For me, simpleITK is the best library for treating or processing Dicom and nifti files because it contains all of the functions required for both formats and even offers options not found in other medical imaging libraries, such as converting nifti files into Dicoms series (this is one of the most sensitive conversions, and not available in other libraries).
Upload a nifti file
Lots of other options are available, you can check them out [here](https://simpleitk.readthedocs.io/en/master/link_DicomSeriesReader_docs.html).

```Python
import SimpleITK as sitk

nifti_file = sitk.ReadImage('The path to the nifti file')
```

---

# VTK

[🔗](https://vtk.org/documentation/) Kitware provides VTK as a library. It includes not only conversion functions but also additional tools for visualization and medical imaging processing! This library is even used by some well-known software, such as 3D Slicer. This means that all of the features available in 3D Slicer can be found and used in VTK, which is fantastic.
Because there are so many things we can do with this library, I'll provide a link to the documentation so you can directly check and find the exact functions/operations required for your projects!

---

# Numpy-stl

[🔗](https://pythonhosted.org/numpy-stl/) This library is not only used for medical imaging, but it plays an important role in medical imaging by assisting in the creation of 3D meshes from segmentation (manual or automatic) or assisting in the conversion of nifti files into STLs, as discussed in this blog post.

![liver_3D](https://user-images.githubusercontent.com/37108394/201542844-e9128fe0-f18d-466d-a086-c89d99b9ce23.gif)

---


# MedPy

[🔗](https://loli.github.io/medpy/) Another SimpleITK-based library for medical imaging that can be used for reading and writing files such as:
Medical formats:

- ITK MetaImage (.mha/.raw, .mhd)
- Neuroimaging Informatics Technology Initiative (NIfTI) (.nia, .nii, .nii.gz, .hdr, .img, .img.gz)
- Analyze (plain, SPM99, SPM2) (.hdr/.img, .img.gz)
- Digital Imaging and Communications in Medicine (DICOM) (.dcm, .dicom)
- Digital Imaging and Communications in Medicine (DICOM) series (<directory>/)
- Nearly Raw Raster Data (Nrrd) (.nrrd, .nhdr)
- Medical Imaging NetCDF (MINC) (.mnc, .MNC)
- Guys Image Processing Lab (GIPL) (.gipl, .gipl.gz)

## Microscopy formats:

- Medical Research Council (MRC) (.mrc, .rec)
- Bio-Rad (.pic, .PIC)
- LSM (Zeiss) microscopy images (.tif, .TIF, .tiff, .TIFF, .lsm, .LSM)
- Stimulate / Signal Data (SDT) (.sdt)


## Visualization formats:

- VTK images (.vtk)


## Other formats:

- Portable Network Graphics (PNG) (.png, .PNG)
- Joint Photographic Experts Group (JPEG) (.jpg, .JPG, .jpeg, .JPEG)
- Tagged Image File Format (TIFF) (.tif, .TIF, .tiff, .TIFF)
- Windows bitmap (.bmp, .BMP)
- Hierarchical Data Format (HDF5) (.h5 , .hdf5 , .he5)
- MSX-DOS Screen-x (.ge4, .ge5)

---

# MONAI

[🔗](https://docs.monai.io/en/stable/) We've arrived at the final one. As they say, "last for best!" The library can be used to implement machine learning algorithms such as image classification and segmentation. MONAI can perform a wide range of operations, from preprocessing data to deploying your final model. I already have some content about it, such as my free 5-hour course on 3D automatic liver segmentation with MONAI.

---

# Full deep learning course

You can join our waitlist [here](https://pycad.co/deep-learning-for-medical-imaging/).

