using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using Heroius.XuAlgrithms.Utility;

namespace Sample
{
    public abstract class SampleItem
    {
        public SampleItem(string name)
        {
            var part = name.Split('.');
            Section = $"{part[0]} : {Mapping.GetSectionDescription(part[0])}";
            Name = $"{part[1]} : {Mapping.GetDescription(name)}";
            Tooltip = Mapping.GetDescription(name);
        }

        public string Section { get; private set; }
        public string Name { get; private set; }
        public string Tooltip { get; private set; }
        
        public abstract string Execute();
    }
}
