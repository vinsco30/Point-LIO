# Point-LIO

**Branch for degenerancy check**
New features added:
1. Item Publish of the complete odometry message (added velocity);
2. Item Added function ``` compute_metrics ``` where Fisher information matrix is computed (it doesn't work with livox lidars);
3. Item Included topic publication of number of tracked points and state-error covariance trace;

## 1. Introduction

**Important notes:**

A. Please make sure the IMU and LiDAR are **Synchronized**, that's important.

B. Please obtain the saturation values of your used IMU (i.e., accelerator and gyroscope), and the units of the accelerator of your used IMU, then modify the .yaml file according to those settings, including values of 'satu_acc', 'satu_gyro', 'acc_norm'. That's improtant.

C. The warning message "Failed to find match for field 'time'." means the timestamps of each LiDAR points are missed in the rosbag file. That is important because Point-LIO processes at the sampling time of each LiDAR point.

D. We recommend to set the **extrinsic_est_en** to false if the extrinsic is given. As for the extrinsic initiallization, please refer to our recent work: [**Robust and Online LiDAR-inertial Initialization**](https://github.com/hku-mars/LiDAR_IMU_Init).

E. If a high odometry output frequency without downsample is required, set ``` publish_odometry_without_downsample ``` as true. Then the warning message of tf "TF_REPEATED_DATA" will pop up in the terminal window, because the time interval between two publish odometery is too small. The following command could be used to suppress this warning to a smaller frequency:

in your catkin_ws/src,

git clone --branch throttle-tf-repeated-data-error git@github.com:BadgerTechnologies/geometry2.git

Then rebuild, source setup.bash, run and then it should be reduced down to once every 10 seconds. If 10 seconds is still too much log output then change the ros::Duration(10.0) to 10000 seconds or whatever you like.

F. If you want to use Point-LIO without imu, set the "imu_en" as false, and provide a predefined value of gavity in "gravity_init" as true as possible in the yaml file, and keep the "use_imu_as_input" as 0.

## **1.1. Developers:**
The codes of this repo are contributed by:
[Dongjiao He (贺东娇)](https://github.com/Joanna-HE) and [Wei Xu (徐威)](https://github.com/XW-HKU)


## **1.2. Related paper**
Our paper is published on Advanced Intelligent Systems(AIS). [Point-LIO](https://onlinelibrary.wiley.com/doi/epdf/10.1002/aisy.202200459), DOI: 10.1002/aisy.202200459

## 2. What can Point-LIO do?
### 2.1 Simultaneous LiDAR localization and mapping (SLAM) without motion distortion

### 2.2 Produce high odometry output frequence and high bandwidth

### 2.3 SLAM with aggressive motions even the IMU is saturated

# **3. Prerequisites**

## **3.1 Ubuntu and [ROS](https://www.ros.org/)**
We tested our code on Ubuntu20.04 with noetic. Ubuntu18.04 and lower versions have problems of environments to support the Point-LIO, try to avoid using Point-LIO in those systems. Additional ROS package is required:
```
sudo apt-get install ros-xxx-pcl-conversions
```

## **3.2 Eigen**
Following the official [Eigen installation](eigen.tuxfamily.org/index.php?title=Main_Page), or directly install Eigen by:
```
sudo apt-get install libeigen3-dev
```

## **3.3 livox_ros_driver**
Follow [livox_ros_driver Installation](https://github.com/Livox-SDK/livox_ros_driver).

*Remarks:*
- Since the Point-LIO supports Livox serials LiDAR, so the **livox_ros_driver** must be installed and **sourced** before run any Point-LIO luanch file.
- How to source? The easiest way is add the line ``` source $Licox_ros_driver_dir$/devel/setup.bash ``` to the end of file ``` ~/.bashrc ```, where ``` $Licox_ros_driver_dir$ ``` is the directory of the livox ros driver workspace (should be the ``` ws_livox ``` directory if you completely followed the livox official document).

## 4. Build
Clone the repository and catkin_make:

```
    cd ~/$A_ROS_DIR$/src
    git clone https://github.com/hku-mars/Point-LIO.git
    cd Point-LIO
    git submodule update --init
    cd ../..
    catkin_make
    source devel/setup.bash
```
- Remember to source the livox_ros_driver before build (follow 3.3 **livox_ros_driver**)
- If you want to use a custom build of PCL, add the following line to ~/.bashrc
```export PCL_ROOT={CUSTOM_PCL_PATH}```

## 5. Directly run

### 5.1 For Avia
Connect to your PC to Livox Avia LiDAR by following  [Livox-ros-driver installation](https://github.com/Livox-SDK/livox_ros_driver), then
```
    cd ~/$Point_LIO_ROS_DIR$
    source devel/setup.bash
    roslaunch point_lio mapping_avia.launch
    roslaunch livox_ros_driver livox_lidar_msg.launch
```
- For livox serials, Point-LIO only support the data collected by the ``` livox_lidar_msg.launch ``` since only its ``` livox_ros_driver/CustomMsg ``` data structure produces the timestamp of each LiDAR point which is very important for Point-LIO. ``` livox_lidar.launch ``` can not produce it right now.
- If you want to change the frame rate, please modify the **publish_freq** parameter in the [livox_lidar_msg.launch](https://github.com/Livox-SDK/livox_ros_driver/blob/master/livox_ros_driver/launch/livox_lidar_msg.launch) of [Livox-ros-driver](https://github.com/Livox-SDK/livox_ros_driver) before make the livox_ros_driver pakage.

### 5.2 For Livox serials with external IMU

mapping_avia.launch theratically supports mid-70, mid-40 or other livox serial LiDAR, but need to setup some parameters befor run:

Edit ``` config/avia.yaml ``` to set the below parameters:

1. LiDAR point cloud topic name: ``` lid_topic ```
2. IMU topic name: ``` imu_topic ```
3. Translational extrinsic: ``` extrinsic_T ```
4. Rotational extrinsic: ``` extrinsic_R ``` (only support rotation matrix)
- The extrinsic parameters in Point-LIO is defined as the LiDAR's pose (position and rotation matrix) in IMU body frame (i.e. the IMU is the base frame). They can be found in the official manual.
5. Saturation value of IMU's accelerator and gyroscope: ```satu_acc```, ```satu_gyro```
6. The norm of IMU's acceleration according to unit of acceleration messages: ``` acc_norm ```

### 5.3 For Velodyne or Ouster (Velodyne as an example)

Step A: Setup before run

Edit ``` config/velodyne.yaml ``` to set the below parameters:

1. LiDAR point cloud topic name: ``` lid_topic ```
2. IMU topic name: ``` imu_topic ``` (both internal and external, 6-aixes or 9-axies are fine)
3. Set the parameter ```timestamp_unit``` based on the unit of **time** (Velodyne) or **t** (Ouster) field in PoindCloud2 rostopic
4. Line number (we tested 16, 32 and 64 line, but not tested 128 or above): ``` scan_line ```
5. Translational extrinsic: ``` extrinsic_T ```
6. Rotational extrinsic: ``` extrinsic_R ``` (only support rotation matrix)
- The extrinsic parameters in Point-LIO is defined as the LiDAR's pose (position and rotation matrix) in IMU body frame (i.e. the IMU is the base frame).
7. Saturation value of IMU's accelerator and gyroscope: ```satu_acc```, ```satu_gyro```
8. The norm of IMU's acceleration according to unit of acceleration messages: ``` acc_norm ```

Step B: Run below
```
    cd ~/$Point_LIO_ROS_DIR$
    source devel/setup.bash
    roslaunch point_lio mapping_velody16.launch
```

Step C: Run LiDAR's ros driver or play rosbag.

### 5.4 PCD file save

Set ``` pcd_save_enable ``` in launchfile to ``` 1 ```. All the scans (in global frame) will be accumulated and saved to the file ``` Point-LIO/PCD/scans.pcd ``` after the Point-LIO is terminated. ```pcl_viewer scans.pcd``` can visualize the point clouds.

*Tips for pcl_viewer:*
- change what to visualize/color by pressing keyboard 1,2,3,4,5 when pcl_viewer is running. 
```
    1 is all random
    2 is X values
    3 is Y values
    4 is Z values
    5 is intensity
```

## 7. Contact us
If you have any questions about this work, please feel free to contact me <hdj65822ATconnect.hku.hk> and Dr. Fu Zhang <fuzhangAThku.hk> via email.
