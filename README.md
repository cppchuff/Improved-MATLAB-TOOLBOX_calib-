＃改进-MATLAB-TOOLBOX_calib-
基于生长的棋盘格角点检测算法改进四点定位过程

主要参考论文：《A Toolbox for Automatic Calibration of Range and Camera Sensors using a single Shot》

MATLAB源码及链接：http://www.cvlibs.net/software/libcbdetect/

TOOLBOX_calib工具箱介绍及源码下载链接：http://www.vision.caltech.edu/bouguetj/calib_doc/

针对原工具箱在针孔相机的棋盘标定后，改进部分算法并加入内参、外参系数，以用于鱼眼相机的标定及畸变矫正。
由于TOOLBOX_calib采用手动四点定位算法，为简化操作，重构相关代码改为基于生长的棋盘格检测方法，改用框选棋盘即可。在预设参数恰当的情况下，能够较为精准的展现棋盘标定情况，并将标定结果传递给后续计算部分。
