// 定义网格大小控制参数 n
n = 10; // 可通过调整此值控制网格大小，n越大网格越密

// 使用 OpenCASCADE 内核进行几何建模，便于使用布尔操作
SetFactory("OpenCASCADE");

// 正方形参数：边长设为 1，中心在原点 (0,0)
square_size = 1;
// 圆形参数：半径为正方形边长的 1/2
circle_radius = square_size / 4; // 直径为边长 1/2，故半径为其 1/4
circle_x = 0; // 圆心 x 坐标
circle_y = 0; // 圆心 y 坐标

// 创建正方形曲面
Rectangle(1) = {-square_size/2, -square_size/2, 0, square_size, square_size};

// 创建圆形曲线（使用完整的圆）
Disk(2) = {circle_x, circle_y, 0, circle_radius, circle_radius};

// 进行布尔差集操作，从正方形中减去圆形区域，形成带孔的正方形
BooleanDifference(3) = { Surface{1}; Delete; } { Surface{2}; Delete; };

// 定义物理曲线（边界条件）
// 获取正方形外边界曲线
Curve Loop(1) = {1, 2, 3, 4}; // 假设布尔操作后外边界曲线标签为 1, 2, 3, 4
Physical Curve("Γ₁") = {1};  // 底部
Physical Curve("Γ₂") = {2};  // 右侧
Physical Curve("Γ₃") = {3};  // 顶部（移动壁面）
Physical Curve("Γ₄") = {4};  // 左侧

// 定义圆形孔洞的边界曲线
// 获取圆形内部边界曲线（布尔操作后生成的新曲线）
Curve Loop(2) = {5}; // 假设圆形边界曲线标签为 5
Physical Curve("Γ₅") = {5}; // 假设圆形边界曲线标签为 5

// 定义计算域的物理曲面（带孔的正方形区域）
Physical Surface("Ω") = {3};

// 设置全局网格尺寸因子
MeshSizeFactor = 1/n; // 根据 n 调整网格大小

// 生成二维网格
Mesh 2;

// 可选：保存网格文件
Save"cylinder_tri_10.msh";