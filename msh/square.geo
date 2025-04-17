// 顶盖驱动方腔几何模型 (Geometry of the lid driven cavity)
// 完全匹配 Figure 2.1 的尺寸和边界条件

// ---------- 几何定义 ----------
// 定义方腔角点（左下角为原点，边长d=0.1 m）
Point(1) = {0,    0,    0, 0.005};  // 左下角 (x=0, y=0)
Point(2) = {0.1,  0,    0, 0.005};  // 右下角 (x=0.1, y=0)
Point(3) = {0.1,  0.1,  0, 0.005};  // 右上角 (x=0.1, y=0.1) → 顶盖
Point(4) = {0,    0.1,  0, 0.005};  // 左上角 (x=0, y=0.1)

// 定义边界线
Line(1) = {1, 2};  // 底部壁面 (Bottom wall)
Line(2) = {2, 3};  // 右侧壁面 (Right wall)
Line(3) = {3, 4};  // 顶部壁面 (Top wall → Lid)
Line(4) = {4, 1};  // 左侧壁面 (Left wall)

// 生成闭合曲线环并创建表面
Curve Loop(1) = {1, 2, 3, 4};  
Plane Surface(1) = {1};

// ---------- 网格控制 ----------  
// 结构化网格划分（四边形单元）
Transfinite Curve {1, 2, 3, 4} = 20;  // 每条边划分20个单元
Transfinite Surface {1} = {1, 2, 3, 4};  
Recombine Surface {1};  // 生成四边形网格

// ---------- 物理组定义 ----------
// 边界条件分组（与OpenFOAM兼容）
Physical Curve("MovingWall") = {3};    // 顶盖 (Ux=1 m/s)
Physical Curve("FixedWalls") = {1, 2, 4};  // 固定壁面
Physical Surface("Fluid") = {1};       // 流体区域