﻿using System;
using System.Collections.Generic;
using System.Linq;

namespace Heroius.XuAlgrithms.Utility
{
    /// <summary>
    /// 提供需求方法和算法名之间的映射
    /// </summary>
    public static class Mapping
    {
        const string ns = "Heroius.XuAlgrithms";

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
            { "TRMUL", "实矩阵相乘" }
        };
    }
}