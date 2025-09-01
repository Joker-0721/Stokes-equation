// 定义网格参数
n = 40;  // 基础划分参数

// 1. 定义几何实体（点→线→面）
Point(1) = {0, 0, 0};
Point(2) = {1/2, 0, 0};
Point(3) = {1/2, 1/2, 0};
Point(4) = {1, 1/2, 0};
Point(5) = {1, 1, 0};
Point(6) = {1/2, 1, 0};
Point(7) = {0, 1, 0};
Point(8) = {0, 1/2, 0};

Line(1) = {1, 2};  // 底部边界
Line(2) = {2, 3};  // 右侧边界
Line(3) = {3, 4};  // 顶部边界
Line(4) = {4, 5};  // 左侧边界S
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {3, 6}; // mid_1
Line(10) = {3, 8}; // mid_2

Curve Loop(1) = {1, 2, 10, 8}; // 左下部分
Curve Loop(2) = {3, 4, 5, -9}; // 右上部分
Curve Loop(3) = {-10, 9, 6, 7}; // 左上部分
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// 2. 配置结构化网格
Transfinite Curve{1,2,3,4,5,6,7,8,9,10} = n/2+1;  // 水平方向9节点
Transfinite Surface{1,2,3};
Recombine Surface{1,2,3};

// 3. 定义物理组（ASCII命名）
Physical Curve("Γ₁") = {1};  // 底部
Physical Curve("Γ₂") = {2};  // 右mid_1侧
Physical Curve("Γ₃") = {3};  // 右mid_2侧
Physical Curve("Γ₄") = {4};  // 右侧
Physical Curve("Γ₅") = {5};  // 顶部
Physical Curve("Γ₆") = {6};  // 左侧
Physical Curve("Γ₇") = {7};  // 左mid_2侧
Physical Curve("Γ₈") = {8};  // 左mid_1侧
Physical Curve("Γ₉") = {9};  // mid_1
Physical Curve("Γ₁₀") = {10}; // mid_2
Physical Surface("Ω") = {1,2,3}; // 计算域

// 4. 网格生成配置
Mesh.Algorithm = 8;          // 使用Frontal算
Mesh.MshFileVersion = 4.1;   // 明确指定格式版本
Geometry.AutoCoherence = 1;  // 自动同步几何实体

// 5. 生成并保存网格
Mesh 2;
Save "test_quad_40.msh";