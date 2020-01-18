using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Windows;
using System.Windows.Controls;
using System.Linq;

namespace Sample
{
    /// <summary>
    /// MainWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
            Samples = new List<SampleSection>();
            foreach (var key in Heroius.XuAlgrithms.Utility.Mapping.GetAllAlgrithmNames())
            {
                try
                {
                    var o = Activator.CreateInstance(null, $"Sample.{key}").Unwrap() as SampleItem;
                    if (!Samples.Exists(s=>s.Section == o.Section))
                    {
                        Samples.Add(new SampleSection() { Section = o.Section, Algrithms = new List<SampleItem>() });
                    }
                    Samples.First(s=>s.Section == o.Section).Algrithms.Add(o);
                }
                catch {
                    //pass undefined sample
                }
            }
            DataContext = this;
        }

        public List<SampleSection> Samples { get; set; }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            var item = (sender as Button).DataContext as SampleItem;
            try
            {
                TbResult.Text = item.Execute();
            }
            catch (Exception ex) { 
                TbResult.Text = ex.Message; }
        }
    }

    public class SampleSection
    {
        public string Section { get; set; }
        public List<SampleItem> Algrithms { get; set; }
    }
}
