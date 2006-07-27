/***************************************************************************
 *            lohner_integrator.tpl
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
//#define DEBUG

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../declarations.h"

#include "../utility/stlio.h"
#include "../base/array.h"
#include "../numeric/arithmetic.h"
#include "../numeric/function.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/interval_vector.h"
#include "../linear_algebra/transformation_system.h"

#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_matrix.h"

#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid.h"
#include "../geometry/grid_set.h"

#include "../system/vector_field.h"
#include "../system/affine_vector_field.h"

#include "../evaluation/integrator.h"

namespace Ariadne {
  namespace Evaluation {
    
    template<typename R>
    class IntervalParallelotope
    {
     public:
      IntervalParallelotope(const LinearAlgebra::IntervalVector<R>& c, const LinearAlgebra::IntervalMatrix<R>& A) : _c(c), _A(A) { }
      
      size_type dimension() const { return _c.size(); }
      
      const LinearAlgebra::IntervalVector<R>& centre() const { return _c; }
      LinearAlgebra::IntervalVector<R>& centre() { return _c; }

      const LinearAlgebra::IntervalMatrix<R>& generators() const { return _A; }
      LinearAlgebra::IntervalMatrix<R>& generators() { return _A; }

      Geometry::Parallelotope<R> over_approximating_parallelotope() const;
      Geometry::Parallelotope<R> over_approximating_parallelotope(R err) const;
     private:
      LinearAlgebra::IntervalVector<R> _c;
      LinearAlgebra::IntervalMatrix<R> _A;
    };
   
    template<typename R>
    inline Geometry::Parallelotope<R> 
    IntervalParallelotope<R>::over_approximating_parallelotope() const 
    {
       return this->over_approximating_parallelotope(R(R(1)/65536));
    }
    
    template<typename R>
    Geometry::Parallelotope<R> 
    IntervalParallelotope<R>::over_approximating_parallelotope(R proposed_err) const 
    {
#ifdef DEBUG
      std::cerr << "IntervalParallelotope<R>::over_approximating_parallelotope() const" << std::endl;
#endif
      typedef typename numerical_traits<R>::field_extension_type F;
      
      size_type n=this->dimension();
      
      LinearAlgebra::IntervalMatrix<R> A=this->generators();
      
      LinearAlgebra::Vector<R> cmid=this->centre().centre();
      LinearAlgebra::Matrix<R> Amid=this->generators().centre();
      
      LinearAlgebra::Matrix<R> D(n,n);
      for(size_type i=0; i!=n; ++i) {
        D(i,i)=this->centre()(i).radius();
      }
      
      LinearAlgebra::Matrix<F> AinvF=LinearAlgebra::inverse(Amid);
      LinearAlgebra::IntervalMatrix<R> Ainv=LinearAlgebra::approximate(AinvF, proposed_err);
      
      R err = upper_norm(Ainv*LinearAlgebra::IntervalMatrix<R>(D+A));
#ifdef DEBUG
      std::cerr << "error=" << err << std::endl;
#endif
      return Geometry::Parallelotope<R>(cmid,err*Amid);
    }
    
    template<typename R>
    std::ostream& 
    operator<<(std::ostream& os, const IntervalParallelotope<R>& ip)
    {
      return os << "IntervalParallelotope(\n  centre=" << ip.centre() << "\n  generators=" << ip.generators() << "\n)\n";
    }
    

    template<typename R>
    class IntervalZonotope
    {
     public:
      IntervalZonotope(const LinearAlgebra::IntervalVector<R>& c, const LinearAlgebra::IntervalMatrix<R>& A) : _c(c), _A(A) { }
      
      size_type dimension() const { return _c.size(); }
      
      const LinearAlgebra::IntervalVector<R>& centre() const { return _c; }
      LinearAlgebra::IntervalVector<R>& centre() { return _c; }

      const LinearAlgebra::IntervalMatrix<R>& generators() const { return _A; }
      LinearAlgebra::IntervalMatrix<R>& generators() { return _A; }

      Geometry::Zonotope<R> over_approximating_zonotope() const;
      Geometry::Zonotope<R> over_approximating_zonotope(R err) const;
     private:
      LinearAlgebra::IntervalVector<R> _c;
      LinearAlgebra::IntervalMatrix<R> _A;
    };
   
    template<typename R>
    inline Geometry::Zonotope<R> 
    IntervalZonotope<R>::over_approximating_zonotope() const 
    {
       return this->over_approximating_zonotope(R(R(1)/65536));
    }
    
    template<typename R>
    Geometry::Zonotope<R> 
    IntervalZonotope<R>::over_approximating_zonotope(R proposed_err) const 
    {
#ifdef DEBUG
      std::cerr << "IntervalZonotope<R>::over_approximating_zonotope() const" << std::endl;
      std::cerr << *this << std::endl;
#endif
      typedef typename numerical_traits<R>::field_extension_type F;
      
      LinearAlgebra::IntervalMatrix<R> A=this->generators();
      
      size_type n=this->dimension();
      size_type m=A.size2();
      
      LinearAlgebra::Vector<R> cmid=this->centre().centre();
      LinearAlgebra::Matrix<R> Amid=this->generators().centre();
      
      LinearAlgebra::Matrix<R> D(n,m);
      for(size_type i=0; i!=std::min(n,m); ++i) {
        D(i,i)=this->centre()(i).radius();
      }
      
      LinearAlgebra::Matrix<F> AinvF=LinearAlgebra::inverse(Amid);
      LinearAlgebra::IntervalMatrix<R> Ainv=LinearAlgebra::approximate(AinvF, proposed_err);
      
      R err = upper_norm(Ainv*LinearAlgebra::IntervalMatrix<R>(D+A));
#ifdef DEBUG
      std::cerr << "error=" << err << std::endl;
#endif
      return Geometry::Zonotope<R>(cmid,err*Amid);
    }
    
    template<typename R>
    std::ostream& 
    operator<<(std::ostream& os, const IntervalZonotope<R>& ip)
    {
      return os << "IntervalZonotope(\n  centre=" << ip.centre() << "\n  generators=" << ip.generators() << "\n)\n";
    }



    template<typename R>
    C1LohnerIntegrator<R>::C1LohnerIntegrator(const R& maximum_step_size, const R& lock_to_grid_time, const R& maximum_basic_set_radius)
      : C1Integrator<R>(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
    {
    }

    template<typename R>
    Geometry::Rectangle<R> 
    C0LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                            const Geometry::Rectangle<R>& initial_set, 
                                            R& step_size) const
    {
#ifdef DEBUG      
      std::cerr << "integrate(const VectorField<R>& vf, const Geometry::Rectangle<R>& r, const R& t)" << std::endl;
#endif

      assert(vector_field.dimension()==initial_set.dimension());
      
      const System::VectorField<R>& vf(vector_field);
      Geometry::Rectangle<R> r=initial_set;
      R& h=step_size;
      
      Geometry::Rectangle<R> q=estimate_flow_bounds(vf,r,h);
      
      LinearAlgebra::IntervalVector<R> fq=vf.apply(q);
      r=r+(h*fq);
#ifdef DEBUG
      std::cerr << "suggested stepsize=" << step_size << std::endl;
                
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << q << std::endl;

      std::cerr << "derivative=" << fq << std::endl;

      std::cerr << "position=" << r << std::endl;
#endif
      return r;
    }

    template<typename R>
    Geometry::Rectangle<R> 
    C0LohnerIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                             const Geometry::Rectangle<R>& initial_set, 
                                             R& step_size) const
    {
#ifdef DEBUG
      std::cerr << "integrate(const VectorField<R>& vf, const Geometry::Rectangle<R>& r, const R& t)" << std::endl;
#endif
      assert(vector_field.dimension()==initial_set.dimension());
      
      const System::VectorField<R>& vf(vector_field);
      Geometry::Rectangle<R> r=initial_set;
      R& h=step_size;

      Geometry::Rectangle<R> q=estimate_flow_bounds(vf,r,h);
      
      LinearAlgebra::IntervalVector<R> fq=vf.apply(q);
      
      r=r+(Interval<R>(R(0),h)*fq);

#ifdef DEBUG
      std::cerr << "suggested stepsize=" << step_size << std::endl;
                
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << q << std::endl;

      std::cerr << "derivative=" << fq << std::endl;
      std::cerr << "position=" << r << std::endl;
#endif

      return r;
    }
    
    template<typename R>
    Geometry::Parallelotope<R> 
    C1LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                            const Geometry::Parallelotope<R>& initial_set, 
                                            R& step_size) const
    {
      const System::VectorField<R>* cvf_ptr=&vector_field;
      const System::AffineVectorField<R>* cavf_ptr=dynamic_cast< const System::AffineVectorField<R>* >(cvf_ptr);

#ifdef DEBUG
      if(cavf_ptr) { std::cerr << typeid(*cavf_ptr).name(); } else { std::cerr << "void"; } std::cerr << std::endl;
      std::cerr << cvf_ptr << "  " << cavf_ptr << "\n" << std::endl;
#endif
      if(cavf_ptr != NULL) {
        return integration_step(*cavf_ptr,initial_set,step_size);
      }

#ifdef DEBUG
      std::cerr << "integration_step(VectorField<R>, Parallelotope<R>, R)\n";
#endif

      typedef typename numerical_traits<R>::field_extension_type F;
       
      assert(vector_field.dimension()==initial_set.dimension());

      const System::VectorField<R>& vf(vector_field);
      Geometry::Parallelotope<R> p=initial_set;
      const size_type n=p.dimension();
      R& h=step_size;
      const LinearAlgebra::Matrix<R> id=LinearAlgebra::identity_matrix<R>(n);
      
      R err=norm(p.generators())/65536;
      
      Geometry::Rectangle<R> b=this->estimate_flow_bounds(vf,p.bounding_box(),h);
#ifdef DEBUG
      std::cerr << "suggested stepsize=" << step_size << std::endl;
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << b << std::endl;
#endif
      
      LinearAlgebra::IntervalVector<R> f=vf.apply(b);
      LinearAlgebra::IntervalMatrix<R> df=vf.derivative(b);
      
      R l=upper_log_norm(df);

#ifdef DEBUG
      std::cerr << "flow=" << f << std::endl;
      std::cerr << "jacobian=" << df << std::endl;
      std::cerr << "jacobian centre=" << df.centre() << std::endl;
      std::cerr << "logarithmic_norm=" << l << std::endl;
#endif

      l=max(l,R(1));
      
      //integration_step(): the varaible l is useless!
      
      LinearAlgebra::IntervalMatrix<R> hdf=h*df;
      LinearAlgebra::IntervalMatrix<R> dphi=exp(hdf);
      
      Geometry::Point<R> c=p.centre();
      Geometry::Rectangle<R> rc=Geometry::Rectangle<R>(c,c);
      Geometry::Rectangle<R> phic=refine_flow_bounds(vf,rc,b,h);
      
      LinearAlgebra::IntervalVector<R> phicv=phic.position_vectors();
      LinearAlgebra::IntervalMatrix<R> zv(dphi*p.generators());

      IntervalParallelotope<R> ip(phicv,zv);
      p=ip.over_approximating_parallelotope();

#ifdef DEBUG
      std::cerr << "stepsize*jacobian=" << hdf << std::endl;
      std::cerr << "flow derivative=" << dphi << std::endl;

      std::cerr << "centre=" << c << std::endl;
      std::cerr << "bounds on centre=" << phic << std::endl;
      std::cerr << "bounds on centre <vector>=" << phicv << std::endl;
      std::cerr << "zonotopic <vector> to add to image of centre=" << zv << std::endl;
      std::cerr << "approximating interval parallelotope=" << ip << std::endl;
      std::cerr << "new approximation=" << p << std::endl;
      std::cerr << "Done integration step\n\n\n" << std::endl;
#endif  
      return p;
    }

    template<typename R>
    Geometry::Parallelotope<R> 
    C1LohnerIntegrator<R>::integration_step(const System::AffineVectorField<R>& vector_field, 
                                            const Geometry::Parallelotope<R>& initial_set, 
                                            R& step_size) const
    {
#ifdef DEBUG
      std::cerr << "integration_step(AffineVectorField<R>, Parallelotope<R>, R)\n";
#endif
      const System::AffineVectorField<R>& vf=vector_field;
      Geometry::Parallelotope<R> p=initial_set;
      R& h=step_size;
      
      R max_error=LinearAlgebra::norm(p.generators())/16;
      assert(max_error>0);
      
#ifdef DEBUG
      std::cerr << "parallelotope generators=" << p.generators() << std::endl;
      std::cerr << "maximum allowed error=" << max_error << std::endl;
      
      std::cerr << "jacobian=" << vf.A() << std::endl;
      std::cerr << "step size=" << h << std::endl;
#endif

      LinearAlgebra::Matrix<R> D=LinearAlgebra::exp_Ah_approx(vf.A(),h,max_error);
      /* Write phi(x)=D x0 + P b */
      LinearAlgebra::Matrix<R> P=LinearAlgebra::exp_Ah_sub_id_div_A_approx(vf.A(),h,max_error);
      
      LinearAlgebra::IntervalMatrix<R> iD(LinearAlgebra::IntervalMatrix<R>(D,max_error));
      LinearAlgebra::IntervalMatrix<R> iP(LinearAlgebra::IntervalMatrix<R>(P,max_error));
      
      IntervalParallelotope<R> img(iD*p.centre().position_vector()+iP*vf.b(),iD*p.generators());
      p=img.over_approximating_parallelotope();      

#ifdef DEBUG
      std::cerr << "twist=" << P << std::endl;
      std::cerr << "approximate derivative=" << D << std::endl;
      
      std::cerr << "approximating derivative=" << iD << std::endl;
      std::cerr << "approximating twist=" << iP << std::endl;
      
      std::cerr << "interval parallelotope=" << img << std::endl;
      std::cerr << "parallelotope=" << p << std::endl;
#endif
      return p;      
    }
    

template<typename R>
    Geometry::Zonotope<R> 
    C1LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                            const Geometry::Zonotope<R>& initial_set, 
                                            R& step_size) const
    {
      const System::VectorField<R>* cvf_ptr=&vector_field;
      const System::AffineVectorField<R>* cavf_ptr=dynamic_cast< const System::AffineVectorField<R>* >(cvf_ptr);

#ifdef DEBUG
      if(cavf_ptr) { std::cerr << typeid(*cavf_ptr).name(); } else { std::cerr << "void"; } std::cerr << std::endl;
      std::cerr << cvf_ptr << "  " << cavf_ptr << "\n" << std::endl;
#endif
      if(cavf_ptr != NULL) {
        return integration_step(*cavf_ptr,initial_set,step_size);
      }

#ifdef DEBUG
      std::cerr << "integration_step(VectorField<R>, Zonotope<R>, R)" << std::endl<<std::flush;
#endif

      typedef typename numerical_traits<R>::field_extension_type F;
       
      using namespace LinearAlgebra;
      using namespace Geometry;
      using namespace System;

      assert(vector_field.dimension()==initial_set.dimension());

      const VectorField<R>& vf(vector_field);
      Zonotope<R> p=initial_set;
      const size_type n=p.dimension();
      R& h=step_size;
      const Matrix<R> id=identity_matrix<R>(n);
      
      R err=norm(p.generators())/65536;
      
      Rectangle<R> b=this->estimate_flow_bounds(vf,p.bounding_box(),h);
#ifdef DEBUG
      std::cerr << "suggested stepsize=" << step_size << std::endl;
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << b << std::endl;
#endif
      
      IntervalVector<R> f=vf.apply(b);
      IntervalMatrix<R> df=vf.derivative(b);
      
      R l=upper_log_norm(df);

#ifdef DEBUG
      std::cerr << "flow=" << f << std::endl;
      std::cerr << "jacobian=" << df << std::endl;
      std::cerr << "jacobian centre=" << df.centre() << std::endl;
      std::cerr << "logarithmic_norm=" << l << std::endl;
#endif

      l=max(l,R(1));
      
      //integration_step(): the varaible l is useless!
      
      IntervalMatrix<R> hdf=h*df;
      IntervalMatrix<R> dphi=exp(hdf);
      
      Point<R> c=p.centre();
      Rectangle<R> rc=Geometry::Rectangle<R>(c,c);
      Rectangle<R> phic=refine_flow_bounds(vf,rc,b,h);
      
      IntervalVector<R> phicv=phic.position_vectors();
      IntervalMatrix<R> zv(dphi*p.generators());

      IntervalZonotope<R> ip(phicv,zv);
      p=ip.over_approximating_zonotope();

#ifdef DEBUG
      std::cerr << "stepsize*jacobian=" << hdf << std::endl;
      std::cerr << "flow derivative=" << dphi << std::endl;

      std::cerr << "centre=" << c << std::endl;
      std::cerr << "bounds on centre=" << phic << std::endl;
      std::cerr << "bounds on centre <vector>=" << phicv << std::endl;
      std::cerr << "zonotopic <vector> to add to image of centre=" << zv << std::endl;
      std::cerr << "approximating interval zonotope=" << ip << std::endl;
      std::cerr << "new approximation=" << p << std::endl;
#endif  

      return p;
    }

    template<typename R>
    Geometry::Zonotope<R> 
    C1LohnerIntegrator<R>::integration_step(const System::AffineVectorField<R>& vector_field, 
                                            const Geometry::Zonotope<R>& initial_set, 
                                            R& step_size) const
    {
#ifdef DEBUG
      std::cerr << "integration_step(AffineVectorField<R>, Parallelotope<R>, R)\n";
#endif
      const System::AffineVectorField<R>& vf=vector_field;
      Geometry::Zonotope<R> p=initial_set;
      R& h=step_size;
      
      R max_error=LinearAlgebra::norm(p.generators())/16;
      assert(max_error>0);
      
#ifdef DEBUG
      std::cerr << "parallelotope generators=" << p.generators() << std::endl;
      std::cerr << "maximum allowed error=" << max_error << std::endl;
      
      std::cerr << "jacobian=" << vf.A() << std::endl;
      std::cerr << "step size=" << h << std::endl;
#endif

      LinearAlgebra::Matrix<R> D=LinearAlgebra::exp_Ah_approx(vf.A(),h,max_error);
      /* Write phi(x)=D x0 + P b */
      LinearAlgebra::Matrix<R> P=LinearAlgebra::exp_Ah_sub_id_div_A_approx(vf.A(),h,max_error);
      
      LinearAlgebra::IntervalMatrix<R> iD(LinearAlgebra::IntervalMatrix<R>(D,max_error));
      LinearAlgebra::IntervalMatrix<R> iP(LinearAlgebra::IntervalMatrix<R>(P,max_error));
      
      IntervalZonotope<R> img(iD*p.centre().position_vector()+iP*vf.b(),iD*p.generators());
      p=img.over_approximating_zonotope();      

#ifdef DEBUG
      std::cerr << "twist=" << P << std::endl;
      std::cerr << "approximate derivative=" << D << std::endl;
      
      std::cerr << "approximating derivative=" << iD << std::endl;
      std::cerr << "approximating twist=" << iP << std::endl;
      
      std::cerr << "interval parallelotope=" << img << std::endl;
      std::cerr << "parallelotope=" << p << std::endl;
#endif
      return p;      
    }

    
    template<typename R>
    Geometry::Zonotope<R> 
    C1LohnerIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                             const Geometry::Zonotope<R>& initial_set, 
                                             R& step_size) const
    {

#ifdef DEBUG
      std::cerr << "reachability_step(VectorField<R>, Zonotope<R>, R)\n";
#endif

      typedef typename numerical_traits<R>::field_extension_type F;

      using namespace LinearAlgebra;
      using namespace Geometry;
      using namespace System;

      assert(vector_field.dimension()==initial_set.dimension());

      const VectorField<R>& vf(vector_field);
      Zonotope<R> z=initial_set;
      const size_type n=z.dimension();
      R h=step_size;
      const Matrix<R> id=identity_matrix<R>(n);
      
      /* Throws exception if we can't find flow bounds for given stepsize. */
      Rectangle<R> b=estimate_flow_bounds(vf,z.bounding_box(),h,256);
      
      IntervalVector<R> f=vf.apply(b);
      IntervalMatrix<R> df=vf.derivative(b);
      
      IntervalMatrix<R> dphi=id+Interval<R>(0,h)*df;
      
      Point<R> c=z.centre();
      Rectangle<R> rc=Geometry::Rectangle<R>(c,c);
      Rectangle<R> phic=refine_flow_bounds(vf,rc,b,R(h/2));
      
      IntervalVector<R> fh=(R(h/2)*f);
      
      TransformationSystem<R> zfh=symmetrise(fh);
      
      Matrix<R> mdf=over_approximation(dphi)*z.generators();
      TransformationSystem<R> zv=zfh+TransformationSystem<R>(Vector<R>(n),mdf);
      
      z=phic+zv;
#ifdef DEBUG
      std::cerr << "suggested stepsize=" << step_size << std::endl;
        
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << b << std::endl;
      
      std::cerr << "flow=" << f << "=" << std::endl;
      std::cerr << "jacobian=" << df << std::endl;
      std::cerr << "flow derivative=" << dphi << std::endl;
        
      std::cerr << "centre=" << c << std::endl;
      std::cerr << "bounds on centre=" << phic << std::endl;
      
      std::cerr << "flow times stepsize=" << fh << std::endl;
      std::cerr << "symmetrised flow=" << zfh << std::endl;
      std::cerr << "over approximating Matrix=" << mdf << std::endl;
      
      std::cerr << "approximating zonotope " << z;
#endif

      return z;
    }

  }
}