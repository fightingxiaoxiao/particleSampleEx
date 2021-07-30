# particleSample

## 简介

用于对OpenFOAM的lagrangian颗粒进行采样。

## 版本适配

* OpenFOAM-v2012

## 安装

创建一个保管第三方代码的文件夹如`extend`，然后从github仓库克隆代码。

```shell
cd $HOME/OpenFOAM/OpenFOAM-<version>
mkdir extend
cd extend
git clone https://github.com/fightingxiaoxiao/particleSample
```

将下载的源码拷贝至应用文件夹（可选，不操作也无所谓）。

```shell
cd $HOME/OpenFOAM/OpenFOAM-<version>/applications/utilities/postProcessing/lagrangian
cp -r $HOME/OpenFOAM/OpenFOAM-<version>/extend/particleSample .
```

进入代码目录，执行`wmake`编译。

```shell
cd particleSample
wmake

```

## 使用

将[particleSampleProperties](particleSampleProperties)放入算例的constant目录。

在算例目录执行`particleSample`。如果是并行算例，执行`mpirun -n <number> particleSample -parallel`。
