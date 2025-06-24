// 定义网格参数
n = 4;  // 基础划分参数

// 1. 定义几何实体（点→线→面）
Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

Line(1) = {1, 2};  // 底部边界
Line(2) = {2, 3};  // 右侧边界
Line(3) = {3, 4};  // 顶部边界
Line(4) = {4, 1};  // 左侧边界

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// 2. 配置结构化网格
Transfinite Curve{1,3} = 2*n+1;  // 水平方向9节点
Transfinite Curve{2,4} = 2*n+1;  // 垂直方向9节点
Transfinite Surface{1};
Recombine Surface{1};

// 3. 定义物理组（ASCII命名）
Physical Curve("Γ₁") = {1};  // 底部
Physical Curve("Γ₂") = {2};  // 右侧
Physical Curve("Γ₃") = {3};  // 顶部（移动壁面）
Physical Curve("Γ₄") = {4};  // 左侧
Physical Surface("Ω") = {1}; // 计算域

// 4. 网格生成配置
Mesh.Algorithm = 8;          // 使用Frontal算法
Mesh.MshFileVersion = 4.1;   // 明确指定格式版本
Geometry.AutoCoherence = 1;  // 自动同步几何实体

// 5. 生成并保存网格
Mesh 2;
Save "mesh.msh";