using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Mathematics.alglib;

namespace Mathematics
{
    /*************************************************************************
    1-мерный интерполянт сплайна
    *************************************************************************/
    public class Spline1dinterpolant : alglibobject
    {
        //
        // Public declarations
        //

        public Spline1dinterpolant()
        {
            _innerobj = new Spline1d.spline1dinterpolant();
        }

        public override alglib.alglibobject make_copy()
        {
            return new Spline1dinterpolant((Spline1d.spline1dinterpolant)_innerobj.make_copy());
        }

        //
        // Although some of declarations below are public, you should not use them
        // They are intended for internal use only
        //
        private Spline1d.spline1dinterpolant _innerobj;
        public Spline1d.spline1dinterpolant innerobj { get { return _innerobj; } }
        public Spline1dinterpolant(Spline1d.spline1dinterpolant obj)
        {
            _innerobj = obj;
        }
    }
}
