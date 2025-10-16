// 定义网格参数
n = 6;  // 基础划分参数

// 1. 定义几何实体（点→线→面）
Point(1) = {0, 0, 0};
Point(2) = {2, 0, 0};
Point(3) = {2, -1, 0};
Point(4) = {10, -1, 0};
Point(5) = {10, 1, 0};
Point(6) = {2, 1, 0};
Point(7) = {0, 1, 0};

// 定义曲线
Line(1) = {1, 2};  
Line(2) = {2, 3};  
Line(3) = {3, 4};  
Line(4) = {4, 5};  
Line(5) = {5, 6};  
Line(6) = {6, 7};
Line(7) = {7, 1};
Line(8) = {2, 6};
// 添加新曲线：从点3到点6
Line(9) = {3, 6};

// 定义曲线循环和表面
Curve Loop(1) = {1, 8, 6, 7};  // 四边形
Curve Loop(2) = {3, 4, 5, -9}; // 现在只有4条曲线：四边形

Plane Surface(1) = {1};
Plane Surface(2) = {2};

// 2. 配置结构化网格
// Transfinite曲线设置：点数量必须为正整数
Transfinite Curve{1,2,6,7,8} = n + 1; // 3 points
Transfinite Curve{3,5} = 4*n + 1;     // 9 points
Transfinite Curve{4} = 2*n + 1;       // 5 points
Transfinite Curve{9} = 2*n + 1;       // 5 points, 与曲线4匹配

// 应用Transfinite和Recombine到两个表面
Transfinite Surface{1};
Transfinite Surface{2};
Recombine Surface{1,2};

// 3. 定义物理组（ASCII命名）
// 根据您的需求调整物理组，确保包含所有边界曲线
Physical Curve("Γ₁") = {7};  // 示例边界组
Physical Curve("Γ₂") = {1,2,3};        // 示例边界组
Physical Curve("Γ₃") = {4};    // 示例边界组
Physical Curve("Γ₄") = {5,6}; // 示例边界组
Physical Surface("Ω") = {1, 2};          // 整个区域

// 4. 网格生成配置
Mesh.Algorithm = 8;          // 使用Frontal算法
Geometry.AutoCoherence = 1;  // 自动同步几何实体

// 5. 生成并保存网格
Mesh 2;
Save "backstep_quad_6.msh";