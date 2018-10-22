using System;
using System.Collections.Generic;
using System.Linq;

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
            return ChsNames.Keys.Search(keyword);
        }
        /// <summary>
        /// 在已筛选的算法类名中进一步筛选
        /// </summary>
        /// <param name="collection">已经筛选的算法类名</param>
        /// <param name="keyword">筛选关键字</param>
        /// <returns></returns>
        public static IEnumerable<string> Search(this IEnumerable<string> collection, string keyword)
        {
            return collection.Where(key => ChsNames[key].Contains(keyword));
        }
        /// <summary>
        /// 获取全部定义的算法名
        /// </summary>
        /// <returns></returns>
        public static IEnumerable<string> GetAllAlgrithmNames()
        {
            return ChsNames.Keys.AsEnumerable();
        }

        /// <summary>
        /// 获取算法描述
        /// </summary>
        /// <param name="key">算法名</param>
        /// <returns></returns>
        public static string GetDescription(string key)
        {
            return ChsNames[key];
        }

        internal static readonly Dictionary<string, string> ChsNames = new Dictionary<string, string>() {
            //矩阵运算
            { "TRMUL", "实矩阵相乘" },
            { "TCMUL", "复矩阵相乘" },
            { "RINV", "一般实矩阵求逆" },
            { "CINV", "一般复矩阵求逆" },
            { "SSGJ", "对称正定矩阵求逆" },
            { "TRCH", "托伯利兹矩阵求逆的特兰持方法" },
            { "SDET", "求一般行列式的值" },
            { "RANK", "求矩阵的秩" },
            { "CHOL", "对称正定矩阵的乔里斯基分解与行列式求值" },
            { "LLUU", "矩阵的三角分解" },
            { "MAQR", "一般实矩阵的QR分解" },
            { "MUAV", "一般实矩阵的奇异值分解" },
            { "GINV", "求广义逆的奇异值分解法" },
            //矩阵特征值与特征向量的计算
            { "SSTQ", "求对称三对角阵的全部特征值与特征向量" },
            { "STRQ", "求实对称矩阵全部特征值与特征向量的豪斯荷尔德变换法" },
            { "HHQR", "求赫申伯格矩阵全部特征值的QR方法" },
            { "BGQR", "求一般实矩阵的全部特征值" },
            { "JCBI", "求实对称矩阵特征值与特征向量的雅可比法" },
            { "JCBJ", "求实对称矩阵特征值与特征向量的雅可比过关法" },
            //线性代数方程组的求解
            { "GAUS", "求解实系数方程组的全选主元高斯消去法" },
            { "", "" },
            { "", "" },
            { "", "" },
            { "", "" },
            { "", "" },
            { "", "" },
            { "", "" },
            { "", "" },
            { "", "" },
            { "", "" },
            { "", "" },
            { "", "" },
        };
    }
}
