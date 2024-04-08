using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;
using FunctionParser;

namespace Метод_деформируемого_многогранника
{
    enum Actions
    {
        Reflection,       //отражение
        Stretching,       //растяжение
        Compression,      //сжатие
        GlobalCompression //глобальное сжатие
    }
    internal class Nelder_Meade_Runner //алгоритм с параметрами
    {
        int n; //размерность
        Expression exp; //функция
        double alpha; //коэффициент отражения
        double beta;  //коэффициент растяжения
        double gamma; //коэффициент сжатия
        double l; //начальное отклонение
        double epsilon; //эпсилон
        Point[] startSimplex;
        public Nelder_Meade_Runner(int n, Expression exp, double epsilon, Point startPoint, double alpha = 1, double beta = 2, double gamma = 0.5, double l = 1) //По начальной точке
        {
            this.n = n;
            this.exp = exp;
            this.epsilon = epsilon;
            this.alpha = alpha;
            this.beta = beta;
            this.gamma = gamma;
            this.l = l;
            GenerateSimplex(startPoint);
        }

        public Nelder_Meade_Runner(int n, Expression exp, double epsilon, Point[] startSimplex, double alpha = 1, double beta = 2, double gamma = 0.5) //По начальному симплексу
        {
            this.n = n;
            this.exp = exp;
            this.epsilon = epsilon;
            this.alpha = alpha;
            this.beta = beta;
            this.gamma = gamma;
            this.startSimplex = Point.Clone(startSimplex);
        }
        void GenerateSimplex(Point startPoint)
        {
            int n = startPoint.dimension;
            Point[] startSimplex = new Point[startPoint.dimension + 1];
            double[] startPointCoords = startPoint.coords;
            startSimplex[0] = new Point(startPointCoords);
            for (int i = 0; i < n; i++)
            {
                startPointCoords[i] += l;
                startSimplex[i+1] = new Point(startPointCoords);
                startPointCoords[i] -= l;                
            }
            this.startSimplex = startSimplex;
        }
        public bool IsReadyToStart() //проверка к старту алгоритма
        {
            //условие 1 (Правильный размер симплекса)
            if (n + 1 != startSimplex.Length) return false;
            //условие 2 (Размерности точек совпадают)
            for (int i = 0; i < startSimplex.Length; i++)
                if ((startSimplex[i] != null) || (n != startSimplex[i].dimension))
                    return false;
            //условие 3 (Корректные параметры алгоритма)
            if ((alpha <= 0) || (beta <= 1) || (gamma <= 0) || (gamma >=1) || (epsilon <= 0)) return false;
            //условие 4 (Начальный симплекс линейно независим)
            if (!IsConvexHull()) return false;
            return true;
        }

        public bool IsConvexHull() //проверка линейной независимости
        {
            Point[] vectors = new Point[n];
            for (int i = 0; i < n; i++)
                vectors[i] = startSimplex[i+1] - startSimplex[0];
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    if (!Point.LinearIndependence(vectors[i], vectors[j])) return false;
            return true;
        }
        public List<Configuration> Run()
        {
            List<Configuration> configurations = new List<Configuration>(100);
            Point.CalculateValue(startSimplex, exp);
            Configuration curConf = CreateConfiguration(startSimplex);
            
            return configurations;
        }
        public Configuration NextConfiguration(Configuration conf)
        {
            Point[] newSimplex = Point.Clone(conf.simplex);
            switch (conf.action)
            {
                case Actions.Reflection:
                    newSimplex[n] = conf.reflectPoint;
                    break;
                case Actions.Stretching:
                    newSimplex[n] = conf.stretchPoint;
                    break;
                case Actions.Compression:
                    newSimplex[n] = conf.сompressPoint;
                    break;
                case Actions.GlobalCompression:
                    GlobalComprerssion(newSimplex, exp);
                    break;
            }
            return CreateConfiguration(newSimplex);
        }
        public Configuration CreateConfiguration(Point[] simplex)
        {
            Point[] newSimplex = Point.Clone(simplex);
            Array.Sort(newSimplex, new ValueIncreasingComparer());
            Point centrePoint = Centre(simplex);
            Point reflectPoint = Reflection(newSimplex, centrePoint, alpha, exp);
            Point stretchPoint = null;
            Point compressPoint = null;
            if ((newSimplex[0].value <= reflectPoint.value) && (reflectPoint.value <= newSimplex[n - 1].value)) //Случай 1 (отражение)
                return new Configuration(newSimplex, Actions.Reflection, centrePoint, reflectPoint, stretchPoint, null);
            else if (reflectPoint.value < newSimplex[0].value) //Случаи 1 и 2 (отражение или растяжение)
            {
                stretchPoint = Stretching(centrePoint, reflectPoint, beta, exp);
                if (stretchPoint.value < reflectPoint.value)
                    return new Configuration(newSimplex, Actions.Stretching, centrePoint, reflectPoint, stretchPoint, compressPoint);
                else
                    return new Configuration(newSimplex, Actions.Reflection, centrePoint, reflectPoint, stretchPoint, compressPoint);
            }
            else//случаи 3 и 4 (сжатие или глобальное сжатие)
            {
                compressPoint = Compression(newSimplex, centrePoint, reflectPoint, gamma, exp);
                if (compressPoint.value < Math.Min(newSimplex[n].value, reflectPoint.value))
                    return new Configuration(newSimplex, Actions.Compression, centrePoint, reflectPoint, stretchPoint, compressPoint);
                else
                    return new Configuration(newSimplex, Actions.GlobalCompression, centrePoint, reflectPoint, stretchPoint, compressPoint);
            }
        }
        bool StopCondition(Configuration conf)
        {
            return !(Dispersion(conf.simplex) < epsilon);
        }
        static double Dispersion(Point[] simplex)
        {
            int n = simplex[0].dimension;
            double disp = 0;
            for (int j = 1; j < n + 1; j++)
                disp += Math.Pow(simplex[j].value - simplex[0].value, 2);
            disp /= n;
            return Math.Sqrt(disp);
        }
        static Point Centre(Point[] simplex)
        {
            int n = simplex[0].dimension;
            double[] point = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    point[i] += simplex[j].coords[i];
                point[i] /= n;
            }
            return new Point(point);
        }
        static Point Reflection(Point[] simplex, Point centre, double alpha, Expression exp)
        {
            int n = simplex[0].dimension;
            double[] point = new double[n];
            for (int i = 0; i < n; i++)
                point[i] = centre.coords[i] + alpha * (centre.coords[i] - simplex[n].coords[i]);
            return new Point(point, exp);
        }
        static Point Stretching(Point centre, Point refPoint, double beta, Expression exp)
        {
            int n = centre.dimension;
            double[] point = new double[n];
            for (int i = 0; i < n; i++)
                point[i] = centre.coords[i] + beta * (refPoint.coords[i] - centre.coords[i]);
            return new Point(point, exp);
        }
        static Point Compression(Point[] simplex, Point centre, Point refPoint, double gamma, Expression exp)
        {
            int n = simplex[0].dimension;
            double[] point = new double[n];
            if (simplex[n].value <= refPoint.value)
                for (int i = 0; i < n; i++)
                    point[i] = centre.coords[i] + gamma * (simplex[n].coords[i] - centre.coords[i]);
            else
                for (int i = 0; i < n; i++)
                    point[i] = centre.coords[i] + gamma * (refPoint.coords[i] - centre.coords[i]);
            return new Point(point, exp);
        }
        static void GlobalComprerssion(Point[] simplex, Expression exp)
        {
            int n = simplex[0].dimension;
            double[] newPoint;
            for (int i = 1; i < n + 1; i++)
            {
                newPoint = new double[n];
                for (int j = 0; j < n; j++)
                    newPoint[j] = 0.5 * (simplex[i].coords[j] + simplex[0].coords[j]);
                simplex[i] = new Point(newPoint, exp);
            }
        }
    }
    internal class Configuration //сюда записываем промежуточные результаты работы
    {
        public readonly Point[] simplex; //симплекс
        public readonly Actions action;  //следущее действие
        public readonly Point centrePoint; //вспомогательные точки
        public readonly Point reflectPoint;
        public readonly Point stretchPoint;
        public readonly Point сompressPoint;

        public Configuration(Point[] simplex,  Actions action, Point centrePoint, Point reflectPoint, Point stretchPoint, Point сompressPoint)
        {
            this.simplex = simplex;
            this.action = action;
            this.centrePoint = centrePoint;
            this.reflectPoint = reflectPoint;
            this.stretchPoint = stretchPoint;
            this.сompressPoint = сompressPoint;
        }
    }
}
