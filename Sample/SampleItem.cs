using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sample
{
    public abstract class SampleItem
    {
        public SampleItem(string name)
        {
            Name = name;
            Tooltip = Heroius.XuAlgrithms.Utility.Mapping.GetDescription(Name);
        }

        public string Name { get; set; }
        public string Tooltip { get; set; }

        public abstract string Execute();
    }
}
