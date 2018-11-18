using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace Heroius.XuAlgrithms.Utility
{
    /// <summary>
    /// 提供需求方法和算法名之间的映射
    /// </summary>
    public static class Mapping
    {
        /// <summary>
        /// 在所有定义的算法中查询功能包含关键字的算法类名
        /// </summary>
        /// <param name="keyword">关键字</param>
        /// <returns></returns>
        public static IEnumerable<string> Search(string keyword)
        {
            return GetAllAlgrithmNames().Search(keyword);
        }
        /// <summary>
        /// 在已筛选的算法类名中进一步筛选
        /// </summary>
        /// <param name="collection">已经筛选的算法类名</param>
        /// <param name="keyword">筛选关键字</param>
        /// <returns></returns>
        public static IEnumerable<string> Search(this IEnumerable<string> collection, string keyword)
        {
            return collection.Where(key => GetDescription(key).Contains(keyword));
        }

        /// <summary>
        /// 获取全部定义的完整限定算法名
        /// </summary>
        /// <returns></returns>
        public static IEnumerable<string> GetAllAlgrithmNames()
        {
            var all = descriptions.Select(tlv1 => tlv1.Item3.Select(tlv2 => $"{tlv1.Item1}.{tlv2.Item1}"));
            var result = all.First();
            foreach (var item in all)
            {
                result = result.Concat(item);
            }
            return result;
        }

        /// <summary>
        /// 获取算法描述
        /// </summary>
        /// <param name="FullName">完整限定的算法名</param>
        /// <returns>算法描述文本</returns>
        public static string GetDescription(string FullName)
        {
            return GetAlgrithmInfo(FullName).Item2;
        }

        /// <summary>
        /// 根据完整限定的算法名获取其信息结构
        /// </summary>
        /// <param name="FullName">完整限定算法名</param>
        /// <returns>包含局部算法名及描述的结构</returns>
        public static Tuple<string, string> GetAlgrithmInfo(string FullName)
        {
            var keyparts = FullName.Split('.');
            return descriptions.First(tlv1 => tlv1.Item1 == keyparts[0]).Item3.First(tlv2 => tlv2.Item1 == keyparts[1]);
        }
        /// <summary>
        /// 根据完整限定的算法名获取其方法反射
        /// </summary>
        /// <param name="FullName"></param>
        /// <returns></returns>
        public static MethodInfo GetAlgrithmReflection(string FullName)
        {
            var keyparts = FullName.Split('.');
            return Type.GetType($"Heroius.XuAlgrithms.{keyparts[0]}").GetMethod(keyparts[1]);
        }

        public static string GetSectionDescription(string SectionName)
        {
            return descriptions.First(t => t.Item1 == SectionName).Item2;
        }

        public static IEnumerable<string> GetSectionAlgrithmNames(string SectionName)
        {
            var sect = descriptions.First(tlv1 => tlv1.Item1 == SectionName);
            return sect.Item3.Select(tlv2 => $"{sect.Item1}.{tlv2.Item1}");
        }

        internal static readonly Tuple<string, string, Tuple<string, string>[]>[] descriptions = new Tuple<string, string, Tuple<string, string>[]>[] {
            new Tuple<string, string, Tuple<string, string>[]>("Matrix", "矩阵运算", new Tuple<string, string>[] {
                new Tuple<string, string>("TRMUL", "实矩阵相乘"),
                new Tuple<string, string>( "TCMUL", "复矩阵相乘" ),
                new Tuple<string, string>( "RINV", "一般实矩阵求逆" ),
                new Tuple<string, string>( "CINV", "一般复矩阵求逆" ),
                new Tuple<string, string>( "SSGJ", "对称正定矩阵求逆" ),
                new Tuple<string, string>( "TRCH", "托伯利兹矩阵求逆的特兰持方法" ),
                new Tuple<string, string>( "SDET", "求一般行列式的值" ),
                new Tuple<string, string>( "RANK", "求矩阵的秩" ),
                new Tuple<string, string>( "CHOL", "对称正定矩阵的乔里斯基分解与行列式求值" ),
                new Tuple<string, string>( "LLUU", "矩阵的三角分解" ),
                new Tuple<string, string>( "MAQR", "一般实矩阵的QR分解" ),
                new Tuple<string, string>( "MUAV", "一般实矩阵的奇异值分解" ),
                new Tuple<string, string>( "GINV", "求广义逆的奇异值分解法" )}),
            new Tuple<string, string, Tuple<string, string>[]>("MatrixEigen", "矩阵特征值与特征向量的计算", new Tuple<string, string>[] {
                new Tuple<string, string>( "SSTQ", "求对称三对角阵的全部特征值与特征向量" ),
                new Tuple<string, string>( "STRQ", "约化对称矩阵为对称三对角阵的豪斯荷尔德变换法" ),
                new Tuple<string, string>( "HHBG", "约化一般实矩阵为赫申伯格矩阵的初等相似变换法" ),
                new Tuple<string, string>( "HHQR", "求赫申伯格矩阵全部特征值的QR方法" ),
                new Tuple<string, string>( "JCBI", "求实对称矩阵特征值与特征向量的雅可比法" ),
                new Tuple<string, string>( "JCBJ", "求实对称矩阵特征值与特征向量的雅可比过关法" )}),
            //todo: BGQR in V4
            new Tuple<string, string, Tuple<string, string>[]>("LinearEquations", "线性代数方程组的求解", new Tuple<string, string>[] {
                new Tuple<string, string>( "GAUS", "求解实系数方程组的全选主元高斯消去法" ),
                new Tuple<string, string>( "GJDN", "求解实系数方程组的全选主元高斯-约当消去法" ),
                new Tuple<string, string>( "CGAS", "求解复系数方程组的全选主元高斯消去法" ),
                new Tuple<string, string>( "CJDN", "求解复系数方程组的全选主元高斯-约当消去法" ),
                new Tuple<string, string>( "TRDE", "求解三对角线方程组的追赶法" ),
                new Tuple<string, string>( "BAND", "求解一般带型方程组" ),
                new Tuple<string, string>( "LDLE", "求解对称方程组的分解法" ),
                new Tuple<string, string>( "CHLK", "求解对称正定方程组的平方根法" ),
                new Tuple<string, string>( "TLVS", "求解托伯利兹方程组的列文逊方法" ),
                new Tuple<string, string>( "GSDL", "高斯-赛德尔迭代法" ),
                new Tuple<string, string>( "GRAD", "求解对称正定方程组的共轭梯度法" ),
                new Tuple<string, string>( "GMQR", "求解线性最小二乘问题的豪斯荷尔德变换法" ),
                new Tuple<string, string>( "GMIV", "求解线性最小二乘问题的广义逆法" ),
                new Tuple<string, string>( "BINT", "求解病态方程组" )})/*,
            new Tuple<string, string, Tuple<string, string>[]>("", "非线性方程与方程组的求解", new Tuple<string, string>[] {
                new Tuple<string, string>( "DHRT", "求非线性方程实根的对分法" ),
                new Tuple<string, string>( "NEWT", "求非线性方程一个实根的牛顿法" ),
                new Tuple<string, string>( "ATKN", "求非线性方程一个实根的埃特金迭代法" ),
                new Tuple<string, string>( "FALSE", "求非线性方程一个实根的试位法" ),
                new Tuple<string, string>( "PQRT", "求非线性方程一个实根的连分式法" ),
                new Tuple<string, string>( "QRRT", "求实系数代数方程全部根的QR方法" ),
                new Tuple<string, string>( "SRRT", "求实系数代数方程全部根的牛顿下山法" ),
                new Tuple<string, string>( "CSRT", "求复系数代数方程全部根的牛顿下山法" ),
                new Tuple<string, string>( "SNSE", "求非线性方程组一组实根的梯度法" ),
                new Tuple<string, string>( "NETN", "求非线性方程组一组实根的拟牛顿法" ),
                new Tuple<string, string>( "NGIN", "求非线性方程组最小二乘解的广义逆法" ),
                new Tuple<string, string>( "MTCL", "求非线性方程一个实根的蒙特卡罗法" ),
                new Tuple<string, string>( "CMTC", "求实函数或复函数方程一个复根的蒙特卡罗法" ),
                new Tuple<string, string>( "NMTC", "求非线性方程组一组实根的蒙特卡罗法" )}),
            new Tuple<string, string, Tuple<string, string>[]>("", "多项式计算", new Tuple<string, string>[] {
                new Tuple<string, string>( "", "一维多项式求值" ),
                new Tuple<string, string>( "", "一维多项式多组求值" ),
                new Tuple<string, string>( "", "二维多项式求值" ),
                new Tuple<string, string>( "", "复系数多项式求值" ),
                new Tuple<string, string>( "", "多项式相乘" ),
                new Tuple<string, string>( "", "复系数多项式相乘" ),
                new Tuple<string, string>( "", "多项式相除" ),
                new Tuple<string, string>( "", "复系数多项式相除" ),
                new Tuple<string, string>( "", "实系数多项式类" ),
                new Tuple<string, string>( "", "复系数多项式类" )}) ,
            new Tuple<string, string, Tuple<string, string>[]>("", "复数运算", new Tuple<string, string>[] {
                new Tuple<string, string>( "", "复数乘法" ),
                new Tuple<string, string>( "", "复数出发" ),
                new Tuple<string, string>( "", "复数乘幂" ),
                new Tuple<string, string>( "", "复数的n次方根" ),
                new Tuple<string, string>( "", "复数指数" ),
                new Tuple<string, string>( "", "复数对数" ),
                new Tuple<string, string>( "", "复数正弦" ),
                new Tuple<string, string>( "", "复数余弦" ),
                new Tuple<string, string>( "", "复数类" )}),
            new Tuple<string, string, Tuple<string, string>[]>("", "随机数的产生", new Tuple<string, string>[] {
                new Tuple<string, string>( "", "产生0~1之间均匀分布的一个随机数" ),
                new Tuple<string, string>( "", "产生0~1之间均匀分布的随机数序列" ),
                new Tuple<string, string>( "", "产生任意区间内均匀分布的一个随机整数" ),
                new Tuple<string, string>( "", "产生任意区间内均匀分布的随机整数序列" ),
                new Tuple<string, string>( "", "产生任意均值与方差的正态分布的一个随机数" ),
                new Tuple<string, string>( "", "产生任意均值与方差的正态分布的随机数序列" )}),
            new Tuple<string, string, Tuple<string, string>[]>("", "插值与逼近", new Tuple<string, string>[] {
                new Tuple<string, string>( "LGR", "一元全区间插值" ),
                new Tuple<string, string>( "LG3", "一元三点插值" ),
                new Tuple<string, string>( "PQS", "连分式插值" ),
                new Tuple<string, string>( "HMT", "埃尔米特插值" ),
                new Tuple<string, string>( "ATK", "埃特金逐步插值" ),
                new Tuple<string, string>( "SPL", "光滑插值" ),
                new Tuple<string, string>( "SPL1", "第一种边界条件的三次样条函数插值、微商与积分" ),
                new Tuple<string, string>( "SPL2", "第二种边界条件的三次样条函数插值、微商与积分" ),
                new Tuple<string, string>( "SPL3", "第三种边界条件的三次样条函数插值、微商与积分" ),
                new Tuple<string, string>( "SLQ3", "二元三点插值" ),
                new Tuple<string, string>( "SLGQ", "二元全区间插值" ),
                new Tuple<string, string>( "PIRL", "最小二乘曲线拟合" ),
                new Tuple<string, string>( "CHIR", "切比雪夫曲线拟合" ),
                new Tuple<string, string>( "REMZ", "最佳一致逼近的里米兹方法" ),
                new Tuple<string, string>( "PIR2", "矩形域的最小二乘曲面拟合" )}),
            new Tuple<string, string, Tuple<string, string>[]>("", "数值积分", new Tuple<string, string>[] {
                new Tuple<string, string>( "FFTS", "变步长梯形求积法" ),
                new Tuple<string, string>( "SIMP", "变步长辛卜生求积法" ),
                new Tuple<string, string>( "FPTS", "自适应梯形求积法" ),
                new Tuple<string, string>( "ROMB", "龙贝格求积法" ),
                new Tuple<string, string>( "FPQG", "计算一维积分的连分式法" ),
                new Tuple<string, string>( "PART", "高振荡函数求积法" ),
                new Tuple<string, string>( "LRGS", "勒让德-高斯求积法" ),
                new Tuple<string, string>( "LAGS", "拉盖尔-高斯求积法" ),
                new Tuple<string, string>( "HMGS", "埃尔米特-高斯求积法" ),
                new Tuple<string, string>( "CBSV", "切比雪夫求积法" ),
                new Tuple<string, string>( "MTCL", "计算一维积分的蒙特卡罗法" ),
                new Tuple<string, string>( "SIM2", "变步长辛卜生二重积分法" ),
                new Tuple<string, string>( "GAUS", "计算多重积分的高斯方法" ),
                new Tuple<string, string>( "PQG2", "计算二重积分的连分式法" ),
                new Tuple<string, string>( "MTML", "计算多重积分的蒙特卡罗法" )}),
            new Tuple<string, string, Tuple<string, string>[]>("", "常微分方程组的求解", new Tuple<string, string>[] {
                new Tuple<string, string>( "", "全区间积分的定步长欧拉方法" ),
                new Tuple<string, string>( "", "积分一步的变步长欧拉方法" ),
                new Tuple<string, string>( "", "全区间积分的维梯方法" ),
                new Tuple<string, string>( "", "全区间积分的定步长龙格-库塔方法" ),
                new Tuple<string, string>( "", "积分一步的变步长龙格-库塔方法" ),
                new Tuple<string, string>( "", "积分一步的变步长基尔方法" ),
                new Tuple<string, string>( "", "全区间积分的变步长默森方法" ),
                new Tuple<string, string>( "", "积分一步的连分式法" ),
                new Tuple<string, string>( "", "全区间积分的双边法" ),
                new Tuple<string, string>( "", "全区间积分的阿当姆斯预报校正法" ),
                new Tuple<string, string>( "", "全区间积分的哈明方法" ),
                new Tuple<string, string>( "", "积分一步的特雷纳方法" ),
                new Tuple<string, string>( "", "积分刚性方程组的吉尔方法" ),
                new Tuple<string, string>( "", "求解二阶微分方程边值问题的差分法" )}),
            new Tuple<string, string, Tuple<string, string>[]>("", "数据处理", new Tuple<string, string>[] {
                new Tuple<string, string>( "", "随机样本分析" ),
                new Tuple<string, string>( "", "一元线性回归分析" ),
                new Tuple<string, string>( "", "多元线性回归分析" ),
                new Tuple<string, string>( "", "逐步回归分析" ),
                new Tuple<string, string>( "", "半对数数据相关" ),
                new Tuple<string, string>( "", "对数数据相关" )}),
            new Tuple<string, string, Tuple<string, string>[]>("", "极值问题的求解", new Tuple<string, string>[] {
                new Tuple<string, string>( "", "一维极值连分式法" ),
                new Tuple<string, string>( "", "n维极值连分式法" ),
                new Tuple<string, string>( "", "不等式约束线性规划问题" ),
                new Tuple<string, string>( "", "求n维极值的单形调优法" ),
                new Tuple<string, string>( "", "求约束条件下n维极值的复形调优法" )}),
            new Tuple<string, string, Tuple<string, string>[]>("", "数学变换与滤波", new Tuple<string, string>[] {
                new Tuple<string, string>( "", "傅里叶级数逼近" ),
                new Tuple<string, string>( "", "快速傅里叶变换" ),
                new Tuple<string, string>( "", "快速沃什变换" ),
                new Tuple<string, string>( "", "五点三次平滑" ),
                new Tuple<string, string>( "", "离散随机线性系统的卡尔曼滤波" ),
                new Tuple<string, string>( "", "α-β-γ滤波" )}),
            new Tuple<string, string, Tuple<string, string>[]>("", "特殊函数的计算", new Tuple<string, string>[] {
                new Tuple<string, string>( "", "伽马函数" ),
                new Tuple<string, string>( "", "不完全伽马函数" ),
                new Tuple<string, string>( "", "误差函数" ),
                new Tuple<string, string>( "", "第一类整数阶贝塞耳函数" ),
                new Tuple<string, string>( "", "第二类整数阶贝塞耳函数" ),
                new Tuple<string, string>( "", "变型第一类整数阶贝塞耳函数" ),
                new Tuple<string, string>( "", "变型第二类整数阶贝塞耳函数" ),
                new Tuple<string, string>( "", "不完全贝塔函数" ),
                new Tuple<string, string>( "", "正态分布函数" ),
                new Tuple<string, string>( "", "t-分布函数" ),
                new Tuple<string, string>( "", "χ²-分布函数" ),
                new Tuple<string, string>( "", "F-分布函数" ),
                new Tuple<string, string>( "", "正弦积分" ),
                new Tuple<string, string>( "", "余弦积分" ),
                new Tuple<string, string>( "", "指数积分" ),
                new Tuple<string, string>( "", "第一类椭圆积分" ),
                new Tuple<string, string>( "", "第二类椭圆积分" )}) ,
            new Tuple<string, string, Tuple<string, string>[]>("", "排序", new Tuple<string, string>[] {
                new Tuple<string, string>( "", "冒泡排序" ),
                new Tuple<string, string>( "", "快速排序" ),
                new Tuple<string, string>( "", "希尔排序" ),
                new Tuple<string, string>( "", "堆排序" ),
                new Tuple<string, string>( "", "结构排序" ),
                new Tuple<string, string>( "", "磁盘文件排序" ),
                new Tuple<string, string>( "", "拓扑分类" )}),
            new Tuple<string, string, Tuple<string, string>[]>("", "查找", new Tuple<string, string>[] {
                new Tuple<string, string>( "", "结构体数组的顺序查找" ),
                new Tuple<string, string>( "", "磁盘随机文本文件的顺序查找" ),
                new Tuple<string, string>( "", "有序数组的对分查找" ),
                new Tuple<string, string>( "", "按关键字成员有序的结构体数组的对分查找" ),
                new Tuple<string, string>( "", "按关键字有序的磁盘随机文本文件的对分查找" ),
                new Tuple<string, string>( "", "磁盘随机文本文件的字符串匹配" )})*/ };
    }
}
