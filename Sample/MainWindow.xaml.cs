using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

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
            Samples = new ObservableCollection<SampleItem>();
            foreach (var key in Heroius.XuAlgrithms.Utility.Mapping.GetAllAlgrithmNames())
            {
                var o = Activator.CreateInstance(null, $"Sample.Samples.{key}").Unwrap() as SampleItem;
                Samples.Add(o);
            }
            DataContext = this;
        }

        public ObservableCollection<SampleItem> Samples { get; set; }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            var item = (sender as Button).DataContext as SampleItem;
            TbResult.Text = item.Execute();
        }
    }
}
