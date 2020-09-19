using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Mathematics.alglib;

namespace Mathematics
{
    public class Lspline : Spline1dCalcs
    {
        /*************************************************************************
        Эта подпрограмма строит интеполянт линейного сплайна

        ВХОДНЫЕ ПАРАМЕТРЫ:
            X   -   узлы сплайна, массив[0..N-1]
            Y   -   значения функци, массив[0..N-1]
            N   -   число точек (опционально):
                    * N>=2
                    * если задано, для построения сплайна используются только первые N точек
                    * если не указано, автоматически определяется по размерам X / Y
                      (length(X) должно быть равно length(Y))

        ВЫХОДНЫЕ ПАРАМЕТРЫ:
            C   -   интерполянт сплайна


        ПОРЯДОК ТОЧЕК:

        Подпрограмма автоматически сортирует точки, поэтому вызывающий может передать несортированный массив.

          -- ALGLIB PROJECT --
             Copyright 24.06.2007 by Bochkanov Sergey
        *************************************************************************/
        public static void spline1dbuildlinear(double[] x, double[] y, int n, out Spline1dinterpolant c)
        {
            c = new Spline1dinterpolant();
            Spline1d.spline1dbuildlinear(x, y, n, c.innerobj);
        }

        public static void spline1dbuildlinear(double[] x, double[] y, out Spline1dinterpolant c)
        {
            int n;
            if ((ap.len(x) != ap.len(y)))
                throw new alglibexception("Ошибка при вызове 'spline1dbuildlinear': похоже, что один из аргументов имеет неправильный размер");
            c = new Spline1dinterpolant();
            n = ap.len(x);
            Spline1d.spline1dbuildlinear(x, y, n, c.innerobj);
        }

    }
}
