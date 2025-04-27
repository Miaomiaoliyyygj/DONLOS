import torch
import numpy as np
import time
import h5py
import cv2
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader
from helper1 import definePsf_torch, resamplingOperator_torch, cnlos_reconstruction5
# 参数定义
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

N = 64*64
width = 1
sett = 0  # 如果N, M, width, bin_resolution改变则设为0
snr = 0.5
z_trim = 0
z_offset = 0
c = 3e8
bin_resolution = 32e-12
sampling_coeff = 1.8


# 加载数据 (需要替换为实际数据加载方式)
data = h5py.File('mask001_data.mat', 'r')
rect_data = np.array(data['out'][:], dtype=np.float32)
rect_data = rect_data.transpose(1, 2, 0)
print(rect_data.shape)  # 假设已转换为npy格式

# 确定分辨率
N = rect_data.shape[0]  # 空间分辨率
M = rect_data.shape[2]  # 时间分辨率
max_range = M * c * bin_resolution  # 直方图的最大范围

# 文件路径
savepath_psf = f'psf_{N}_{M}_w{width}.pt'
savepath_mtx = f'mtx_{M}.pt'

if sett == 0:
    # 需要实现definePsf和resamplingOperator的PyTorch版本
    psf = definePsf_torch(N, M, width / max_range)
    mtx, mtxi = resamplingOperator_torch(M)

    # 转换所有 ndarray 到 torch.Tensor
    rect_data_tensor = torch.from_numpy(rect_data).float()  # 确保是 float32
    mtx_tensor = torch.from_numpy(mtx).float()
    mtxi_tensor = torch.from_numpy(mtxi).float()
    psf_tensor = torch.from_numpy(psf).float()

    # 调用函数
    vol, tic_x, tic_y, tic_z = cnlos_reconstruction5(
        rect_data_tensor, width, z_trim, z_offset,
        psf_tensor, mtx_tensor, mtxi_tensor, bin_resolution
    )

    # 调整方向
    vol = torch.flip(vol, [0])
    vol = torch.flip(vol, [2])
    volumn_MxNxN = vol.numpy()
    # zdim = volumn_MxNxN.shape[0] * 100 // 128
    # volumn_MxNxN = volumn_MxNxN[:zdim]
    print('volumn min, %f' % volumn_MxNxN.min())
    print('volumn max, %f' % volumn_MxNxN.max())

    volumn_MxNxN[volumn_MxNxN < 0] = 0
    front_view = np.max(volumn_MxNxN, axis=0)
    print(front_view.shape)
    front_view_new = (front_view / np.max(front_view)) * 255
    front_view_new = cv2.resize(front_view_new, (256, 256))
    cv2.imwrite("re_rsd_1.png", front_view_new)

    y_view = np.max(volumn_MxNxN, axis=1)
    y_view = (y_view / np.max(y_view)) * 255
    y_view = cv2.resize(y_view, (256, 256))
    cv2.imwrite("re_rsd_2.png", y_view)

    x_view = np.max(volumn_MxNxN, axis=2)
    x_view = (x_view / np.max(x_view)) * 255
    x_view = cv2.resize(x_view, (256, 256))
    cv2.imwrite("re_rsd_3.png", x_view)

    volumn_ZxYxX = volumn_MxNxN
    volumn_ZxYxX = volumn_ZxYxX / np.max(volumn_ZxYxX)



