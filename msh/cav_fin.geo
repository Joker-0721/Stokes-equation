// 定义网格参数
n = 4;  // 基础划分参数

// 1. 定义四边形角点
Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

// 2. 定义边界线（确保闭合）
Line(1) = {1, 2};  // 底部
Line(2) = {2, 3};  // 右侧
Line(3) = {3, 4};  // 顶部（移动壁面）
Line(4) = {4, 1};  // 左侧

// 3. 生成闭合表面
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// 4. 设置网格尺寸参数
Mesh.CharacteristicLengthFactor = 1.0/n;  // 控制全局网格密度
Mesh.CharacteristicLengthMin = 0.1/n;
Mesh.CharacteristicLengthMax = 0.5/n;

// 5. 物理组定义
Physical Curve("Γ₁") = {1};  // 底部
Physical Curve("Γ₂") = {2};  // 右侧
Physical Curve("Γ₃") = {3};  // 顶部（移动壁面）
Physical Curve("Γ₄") = {4};  // 左侧
Physical Surface("Ω") = {1}; // 计算域

// 6. 网格生成设置
Mesh.Algorithm = 6;          // 使用前向算法生成三角形网格
Mesh.MshFileVersion = 4.1;   // 设置网格文件版本
Mesh.ElementOrder = 1;        // 线性单元（一阶单元）
Mesh 2;                      // 生成二维三角形网格

// 7. 保存所有数据
Mesh.SaveAll = 1;