# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 15:05:46 2022

@author: trmabdelaal
"""

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import imageio as iio
from tqdm import tqdm

script_dir = os.path.dirname(os.path.abspath(__file__))
print(script_dir)

#### Define functions

def Binarize255(data):
    return (data > 1).astype(np.int_)*255

def Binarize(data):
    return (data > 1).astype(np.int_)

def Highlight_one_cell(data, cell_num):
    return (data == cell_num).astype(np.int_)*255

def rotate_180(data):
    return np.flip(np.flip(data,axis=0),axis=1)

def clip_image(data,x_neg_offset,x_pos_offset,y_neg_offset,y_pos_offset):
    return data[y_neg_offset:(data.shape[0]-y_pos_offset),x_neg_offset:(data.shape[1]-x_pos_offset)]

def High_to_low(pixel_x,pixel_y):
    return int(pixel_x / 5), int(pixel_y / 5)

### Set variables ###
Image_name = 'Tissue_A_02'
IMC_clip_params = (0, 60, 0, 295)
MSI_clip_params = (9, 0, 53, 0)



IMC_images_path = os.path.join(script_dir, "IMC_images+masks", Image_name, "channels")
IMC_masks_path = os.path.join(script_dir, "IMC_images+masks", Image_name)
MSI_images_path = os.path.join(script_dir, "MALDI_singleTIFF", Image_name)

# IMC data (test with one sample "Vimentin")
IMC_img = iio.imread(os.path.join(IMC_images_path, '194Pt_Vimentin.ome.tiff'))


plt.figure()
plt.imshow(IMC_img,cmap="gray",vmax=1)
plt.title("IMC high resolution")
plt.axis('off')
plt.savefig("IMC high resolution.tiff", bbox_inches='tight', dpi=150)
plt.show()

IMC_img = clip_image(IMC_img, *IMC_clip_params)

plt.figure()
plt.imshow(IMC_img,cmap="gray",vmax=1)
plt.title("IMC high resolution clipped")
plt.axis('off')
plt.savefig("IMC high resolution clipped.tiff", bbox_inches='tight', dpi=150)
plt.show()

# LR IMC
IMC_LR = np.zeros((int(IMC_img.shape[0]/5),int(IMC_img.shape[1]/5)))
#IMC_img_binary = Binarize(IMC_img)
for i in range(IMC_LR.shape[0]):
    for j in range(IMC_LR.shape[1]):
        # IMC_LR[i,j] = np.mean(IMC_img_binary[(i*5):(i*5)+5,(j*5):(j*5)+5])
        IMC_LR[i,j] = np.mean(IMC_img[(i*5):(i*5)+5,(j*5):(j*5)+5])
#del i,j,IMC_img_binary
del i,j

plt.figure()
plt.imshow(IMC_LR,cmap="gray",vmax=1)
plt.title("IMC low resolution")
plt.axis('off')
plt.savefig("IMC low resolution.tiff", bbox_inches='tight', dpi=150)
plt.show()

# cell ids
cell_masks = iio.imread(os.path.join(IMC_masks_path, Image_name + '.ome_mask_1.tiff'))
cell_masks = clip_image(cell_masks, *IMC_clip_params)

plt.figure()
plt.imshow(cell_masks,cmap="gray")
plt.title("Cell masks IMC")
plt.axis('off')
plt.show()

plt.figure()
plt.imshow(Binarize255(cell_masks),cmap="gray")
plt.title("Binarized cell masks IMC")
plt.axis('off')
plt.show()

cell_masking = np.array(cell_masks)
unique_cell_idx, counts = np.unique(cell_masking.flatten(), return_counts=True)
del cell_masking

# cell_num=4274
# Highlight_mask = Highlight_one_cell(cell_masks,cell_num)

# plt.figure()
# plt.imshow(Highlight_mask,cmap="gray")
# plt.axis('off')
# plt.show()

# High resolution IMC data
image_list = os.listdir(IMC_images_path)
Channel_name = image_list
for i in range(len(image_list)):
    Channel_name[i] = image_list[i].split('.')[0]
del i,image_list

Pixel_data_HR = pd.DataFrame(0,columns=np.concatenate((['X','Y','Cell_idx'],Channel_name)),
                             index=['{},{}'.format(i,j) for i in range(cell_masks.shape[0]) for j in range(cell_masks.shape[1])])
Pixel_data_HR['X'] = [i for i in range(cell_masks.shape[0]) for _ in range(cell_masks.shape[1])]
Pixel_data_HR['Y'] = [j for _ in range(cell_masks.shape[0]) for j in range(cell_masks.shape[1])]
Pixel_data_HR['Cell_idx'] = [cell_masks[i,j] for i in range(cell_masks.shape[0]) for j in range(cell_masks.shape[1])]
for c in tqdm(Channel_name):
    data = iio.imread(os.path.join('{}/{}.ome.tiff'.format(IMC_images_path,c)))
    data = clip_image(data, *IMC_clip_params)
    Pixel_data_HR[c] = [data[i,j] for i in range(cell_masks.shape[0]) for j in range(cell_masks.shape[1])]
del c, Channel_name
    
IMC_Cell_Expr = Pixel_data_HR.groupby("Cell_idx").mean()
IMC_Cell_Expr = IMC_Cell_Expr.loc[1:,:] # remove cell id=0 (background)

fig = plt.figure()
ax = fig.add_subplot()
plt.scatter(IMC_Cell_Expr["Y"],IMC_Cell_Expr["X"],c=IMC_Cell_Expr["194Pt_Vimentin"],cmap="gray",s=2,vmax=1)
plt.gca().invert_yaxis()
#plt.axis('off')
ax.set_aspect('equal')
ax.set_facecolor("black")
plt.savefig("IMC cell image.tiff", bbox_inches='tight', dpi=150)
plt.show()
del fig,ax

test_image = iio.imread(os.path.join(IMC_images_path, '194Pt_Vimentin.ome.tiff'))

plt.figure()
plt.imshow(clip_image(test_image, *IMC_clip_params),cmap="gray",vmax=1)
plt.title("test IMC image")
plt.axis('off')
plt.savefig("IMC original image.tiff", bbox_inches='tight', dpi=150)
plt.show()

Pixel_data_HR.to_csv("IMC_Pixel_data.csv")
IMC_Cell_Expr.to_csv("IMC_cell_data.csv")

# MSI data (test with one sample "tissue_A_01_0058")
MSI_img = iio.imread(os.path.join(MSI_images_path, Image_name + '_0058.tif'))


plt.figure()
plt.imshow(MSI_img,cmap="gray")
plt.title("MSI")
plt.axis('off')
plt.savefig("MSI.tiff", bbox_inches='tight', dpi=150)
plt.show()

plt.figure()
plt.imshow(rotate_180(MSI_img),cmap="gray")
plt.title("MSI rotated")
plt.axis('off')
plt.savefig("MSI rotated.tiff", bbox_inches='tight', dpi=150)
plt.show()

plt.figure()
plt.imshow(clip_image(rotate_180(MSI_img), *MSI_clip_params),cmap="gray")
plt.title("MSI rotated clipped")
plt.axis('off')
plt.savefig("MSI rotated clipped.tiff", bbox_inches='tight', dpi=150)
plt.show()

## High_to_low mapping
# for i in range(10):
#     for j in range(10):
#         x,y = High_to_low(i,j)
#         print("X = {} and Y = {}".format(x,y))
        
# Cell_idx_LR = pd.DataFrame(columns=['X','Y','Cell_idx'], index=['{},{}'.format(i,j) for i in range(IMC_LR.shape[0]) for j in range(IMC_LR.shape[1])])
# Cell_idx_LR["X"] = [i for i in range(IMC_LR.shape[0]) for _ in range(IMC_LR.shape[1])]
# Cell_idx_LR["Y"] = [j for _ in range(IMC_LR.shape[0]) for j in range(IMC_LR.shape[1])]
# Cell_idx_LR["Cell_idx"] = [list([]) for _ in range(IMC_LR.shape[0]) for _ in range(IMC_LR.shape[1])]

## High_to_low MSI mapping
image_list = os.listdir(MSI_images_path)
Channel_name = image_list
for i in range(len(image_list)):
    Channel_name[i] = image_list[i].split('_')[3][:-4]
del i,image_list

Pixel_data_LR = pd.DataFrame(columns=np.concatenate((['X','Y','Cell_idx'],Channel_name)), index=['{},{}'.format(i,j) for i in range(IMC_LR.shape[0]) for j in range(IMC_LR.shape[1])])
for i in tqdm(range(IMC_LR.shape[0])):
    for j in range(IMC_LR.shape[1]):
        Pixel_data_LR.loc['{},{}'.format(i,j),'X']=i
        Pixel_data_LR.loc['{},{}'.format(i,j),'Y']=j
        Pixel_data_LR.loc['{},{}'.format(i,j),'Cell_idx']=[]
        
for i in tqdm(range(cell_masks.shape[0])):
    for j in range(cell_masks.shape[1]):
        x,y = High_to_low(i, j)
        temp_list = Pixel_data_LR.loc['{},{}'.format(x,y),'Cell_idx']
        if (cell_masks[i,j] != 0):
            temp_list.append(cell_masks[i,j])
        Pixel_data_LR.loc['{},{}'.format(x,y),'Cell_idx'] = temp_list
        
Overlapping_cell_stats=Pixel_data_LR["Cell_idx"].apply(lambda x: len(np.unique(x))).value_counts()

plt.figure()
plt.bar(Overlapping_cell_stats.index,Overlapping_cell_stats)
plt.title("LR pixel stats")
plt.xlabel("#cells mapping to one pixel")
plt.ylabel("count")
plt.savefig("Low resolution stats1.tiff", bbox_inches='tight', dpi=150)
plt.show()

for c in tqdm(Channel_name):
    data = iio.imread(os.path.join(MSI_images_path, Image_name +'_{}.tif'.format(c)))  
    data = clip_image(rotate_180(data), *MSI_clip_params)
    Pixel_data_LR[c] = [data[i,j] for i in range(IMC_LR.shape[0]) for j in range(IMC_LR.shape[1])]
del c, Channel_name
        
LR_pixel_stats = Pixel_data_LR["Cell_idx"].str.len().value_counts()

plt.figure()
plt.bar(LR_pixel_stats.index,LR_pixel_stats)
plt.title("LR pixel stats2")
plt.xlabel("#non-background pixels mapping to one pixel")
plt.ylabel("count")
plt.savefig("Low resolution stats2.tiff", bbox_inches='tight', dpi=150)
plt.show()

MSI_Cell_Expr = Pixel_data_LR.loc[Pixel_data_LR["Cell_idx"].str.len() != 0,:].explode("Cell_idx").groupby("Cell_idx").mean(numeric_only=False)
MSI_Cell_Expr_pure_pixels = Pixel_data_LR.loc[Pixel_data_LR["Cell_idx"].apply(lambda x: len(np.unique(x))) == 1,:].explode("Cell_idx").groupby("Cell_idx").mean(numeric_only=False)

fig = plt.figure()
ax = fig.add_subplot()
color_scale = np.array(MSI_Cell_Expr["0058"])
#color_scale[color_scale>1] = 1
plt.scatter(MSI_Cell_Expr["Y"],MSI_Cell_Expr["X"],c=color_scale,cmap="gray", s=2)
plt.gca().invert_yaxis()
#plt.axis('off')
ax.set_aspect('equal')
ax.set_facecolor("black")
plt.savefig("MSI cell image.tiff", bbox_inches='tight', dpi=150)
plt.show()
del fig,ax

fig = plt.figure()
ax = fig.add_subplot()
color_scale = np.array(MSI_Cell_Expr_pure_pixels["0058"])
#color_scale[color_scale>1] = 1
plt.scatter(MSI_Cell_Expr_pure_pixels["Y"],MSI_Cell_Expr_pure_pixels["X"],c=color_scale,cmap="gray", s=2)
plt.gca().invert_yaxis()
#plt.axis('off')
ax.set_aspect('equal')
ax.set_facecolor("black")
plt.savefig("MSI pure cell image.tiff", bbox_inches='tight', dpi=150)
plt.show()
del fig,ax

test_image = iio.imread(os.path.join(MSI_images_path, Image_name + '_0058.tif')) 
plt.figure()
plt.imshow(clip_image(rotate_180(test_image), *MSI_clip_params),cmap="gray")
plt.title("MSI rotated clipped")
plt.axis('off')
plt.savefig("MSI original image.tiff", bbox_inches='tight', dpi=150)
plt.show()

Pixel_data_HR.to_csv(os.path.join(script_dir, Image_name + "_IMC_Pixel_data.csv"))
MSI_Cell_Expr.to_csv(os.path.join(script_dir, Image_name + "_MSI_cell_data.csv"))
MSI_Cell_Expr_pure_pixels.to_csv(os.path.join(script_dir, Image_name + "_MSI_pure_cell_data.csv"))