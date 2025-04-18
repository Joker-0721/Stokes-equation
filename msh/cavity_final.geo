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

// 4. 结构化网格命令
Transfinite Curve{1,3} = 2*n+1;  // 水平方向9个节点
Transfinite Curve{2,4} = 2*n+1;  // 垂直方向9个节点
Transfinite Surface{1};
Recombine Surface {1};

// 5. 物理组定义
Physical Curve("Γᵗ") = {3};        // 移动壁面（受力边界）
Physical Curve("Γᵍ") = {1, 2, 4};  // 固定壁面（几何约束边界）
Physical Surface("Ω") = {1};        // 计算域

// 6. 网格生成设置
Mesh.Algorithm = 8;     // 选择网格生成算法
Mesh.MshFileVersion = 2;    // 设置网格文件版本
Mesh 2;                 // 生成二维网格
RecombineMesh;          // 重组网格为四边形

// 7. 保存所有数据
Mesh.SaveAll = 1;