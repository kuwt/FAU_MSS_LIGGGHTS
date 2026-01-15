/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef NORMAL_MODEL
NORMAL_MODEL(HERTZ_STIFFNESSLUBRICANT,hertz/stiffnesslubricant,15)
#else
#ifndef NORMAL_MODEL_HERTZ_STIFFNESSLUBRICANT_H_
#define NORMAL_MODEL_HERTZ_STIFFNESSLUBRICANT_H_
#include "contact_models.h"
#include "normal_model_base.h"
#include "global_properties.h"
#include <cmath>

namespace MODEL_PARAMS
{
    inline static ScalarProperty* createMaxSeparationDistanceRatioLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* maxSeparationDistanceRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "maxSeparationDistanceRatio", caller);
      return maxSeparationDistanceRatioScalar;
    }

    inline static ScalarProperty* createFluidViscosityLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* fluidViscosityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidviscosity", caller);
      return fluidViscosityScalar;
    }
}


namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<HERTZ_STIFFNESSLUBRICANT> : public NormalModelBase
  {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
      NormalModelBase(lmp, hsetup, c),
      k_n(NULL),
      k_t(NULL),
      gamma_n(NULL),
      gamma_t(NULL),
      tangential_damping(false),
      limitForce(false),
      displayedSettings(false)
    {
      
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("k_n", &MODEL_PARAMS::createKn);
      registry.registerProperty("k_t", &MODEL_PARAMS::createKt);
      registry.registerProperty("gamma_n", &MODEL_PARAMS::createGamman);
      registry.registerProperty("gamma_t", &MODEL_PARAMS::createGammat);
      registry.registerProperty("maxSeparationDistanceRatio", &MODEL_PARAMS::createMaxSeparationDistanceRatioLub);
      registry.registerProperty("fluidviscosity", &MODEL_PARAMS::createFluidViscosityLub);

    
      registry.connect("k_n", k_n,"model hertz/stiffnesslubricant");
      registry.connect("k_t", k_t,"model hertz/stiffnesslubricant");
      registry.connect("gamma_n", gamma_n,"model hertz/stiffnesslubricant");
      registry.connect("gamma_t", gamma_t,"model hertz/stiffnesslubricant");
      double maxSeparationDistanceRatio = 1.0;
      registry.connect("maxSeparationDistanceRatio", maxSeparationDistanceRatio,"model hertz/stiffnesslubricant");
      registry.connect("fluidviscosity", fluidDynamicViscosity,"model hertz/stiffnesslubricant");
        
      // error checks on coarsegraining
      if(force->cg_active()) {
        error->cg(FLERR,"model hertz/stiffness");
      }

      neighbor->register_contact_dist_factor(maxSeparationDistanceRatio*1.1); 
      if(maxSeparationDistanceRatio < 1.0) {
            error->one(FLERR,"\n\ncohesion model  requires maxSeparationDistanceRatio >= 1.0. Please increase this value.\n");
      }
    }

    // effective exponent for stress-strain relationship
    
    inline double stressStrainExponent()
    {
      return 1.5;
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double meff = sidata.meff;
      double reff = sidata.is_wall ? sidata.radi : (sidata.radi*sidata.radj/(sidata.radi+sidata.radj));
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical && atom->superquadric_flag) {
          reff = sidata.reff;
      }
#endif

      const double polyhertz = sqrt(reff*sidata.deltan);
      double kn = polyhertz*k_n[itype][jtype];
      double kt = polyhertz*k_t[itype][jtype];
      const double gamman = polyhertz*meff*gamma_n[itype][jtype];
      const double gammat = tangential_damping ? polyhertz*meff*gamma_t[itype][jtype] : 0.0;

      if(!displayedSettings)
      {
        displayedSettings = true;

        /*
        if(limitForce)
            if(0 == comm->me) fprintf(screen," NormalModel<HERTZ_STIFFNESS>: will limit normal force.\n");
        */
      }
      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double Fn_damping = -gamman*sidata.vn;
      const double Fn_contact = kn*sidata.deltan;
      double Fn               = Fn_damping + Fn_contact;

      //limit force to avoid the artefact of negative repulsion force
      if(limitForce && (Fn<0.0) )
      {
          Fn = 0.0;
      }

      sidata.Fn = Fn;
      sidata.kn = kn;
      sidata.kt = kt;
      sidata.gamman = gamman;
      sidata.gammat = gammat;

      #ifdef NONSPHERICAL_ACTIVE_FLAG
          double torque_i[3] = {0.0, 0.0, 0.0}; //initialized here with zeros to avoid compiler warnings
          double Fn_i[3] = { Fn * sidata.en[0], Fn * sidata.en[1], Fn * sidata.en[2]};
          if(sidata.is_non_spherical) {
            double xci[3];
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.i], xci);
            vectorCross3D(xci, Fn_i, torque_i);
          }
      #endif
      // apply normal force
      if(sidata.is_wall) {
        const double Fn_ = Fn * sidata.area_ratio;
        i_forces.delta_F[0] += Fn_ * sidata.en[0];
        i_forces.delta_F[1] += Fn_ * sidata.en[1];
        i_forces.delta_F[2] += Fn_ * sidata.en[2];
        #ifdef NONSPHERICAL_ACTIVE_FLAG
                if(sidata.is_non_spherical) {
                  //for non-spherical particles normal force can produce torque!
                  i_forces.delta_torque[0] += torque_i[0];
                  i_forces.delta_torque[1] += torque_i[1];
                  i_forces.delta_torque[2] += torque_i[2];
                }
        #endif
      } else {
        i_forces.delta_F[0] += sidata.Fn * sidata.en[0];
        i_forces.delta_F[1] += sidata.Fn * sidata.en[1];
        i_forces.delta_F[2] += sidata.Fn * sidata.en[2];

        j_forces.delta_F[0] += -i_forces.delta_F[0];
        j_forces.delta_F[1] += -i_forces.delta_F[1];
        j_forces.delta_F[2] += -i_forces.delta_F[2];
        #ifdef NONSPHERICAL_ACTIVE_FLAG
                if(sidata.is_non_spherical) {
                  //for non-spherical particles normal force can produce torque!
                  double xcj[3], torque_j[3];
                  double Fn_j[3] = { -Fn_i[0], -Fn_i[1], -Fn_i[2]};
                  vectorSubtract3D(sidata.contact_point, atom->x[sidata.j], xcj);
                  vectorCross3D(xcj, Fn_j, torque_j);

                  i_forces.delta_torque[0] += torque_i[0];
                  i_forces.delta_torque[1] += torque_i[1];
                  i_forces.delta_torque[2] += torque_i[2];

                  j_forces.delta_torque[0] += torque_j[0];
                  j_forces.delta_torque[1] += torque_j[1];
                  j_forces.delta_torque[2] += torque_j[2];
                }
        #endif
      }
    }

    void computeLubricationForce(
        double ri, double rj,
        const double xi[3], const double xj[3],
        const double vi[3], const double vj[3],
        const double wi[3], const double wj[3],
        double eta,     // fluid viscosity
        double Epsi,    //2h/(r_i,r_j)
        double beta,
        double F_lub[3],   // output: force on particle i
        double t_lub[3]   // output: torque on particle i
    ) {
        double pi=3.1415926;
        // Lubrication coefficients
        double X11A = 6.0*M_PI*ri*(2.0*pow(beta,2)/(pow(1.0+beta,3)*Epsi)+beta*(1.0+7.0*beta+pow(beta,2))/(5.0*pow(1.0+beta,3))*log(1.0/Epsi));
        double Y11A = 6.0*M_PI*ri*(4.0*beta*(1.0+7.0*beta+pow(beta,2))/(15.0*pow(1.0+beta,3)))*log(1.0/Epsi);
        double Y11B = -4.0*M_PI*pow(ri,2)*(beta*(4.0+beta)/(5.0*pow(1.0+beta,2)))*log(1.0/Epsi);
        double Y21B = -4.0*M_PI*pow(rj,2)*((1.0/beta)*(4.0+1.0/beta)/(5.0*pow(1.0+1.0/beta,2)))*log(1.0/Epsi);
        double Y11C = 8.0*M_PI*pow(ri,3)*(2.0*beta/(5.0*(1.0+beta)))*log(1.0/Epsi);
        double Y12C = 8.0*M_PI*pow(ri,3)*(beta/(10.0*(1.0+beta)))*log(1.0/Epsi);

        // --------------------------------
        // Compute n_ij (from j to i)
        // --------------------------------
        double rij[3] = {
            xi[0] - xj[0],
            xi[1] - xj[1],
            xi[2] - xj[2]
        };

        double r = std::sqrt(
            rij[0]*rij[0] +
            rij[1]*rij[1] +
            rij[2]*rij[2]
        );

        double nij[3] = {
            rij[0] / r,
            rij[1] / r,
            rij[2] / r
        };

        // --------------------------------
        // Relative velocity (from j to i)
        // --------------------------------
        double vij[3] = {
            vi[0] - vj[0],
            vi[1] - vj[1],
            vi[2] - vj[2]
        };

        // --------------------------------
        // Normal / tangential velocity split
        // --------------------------------
        double v_n_mag =vij[0]*nij[0] + vij[1]*nij[1] + vij[2]*nij[2];

        double v_n[3] = {
            v_n_mag * nij[0],
            v_n_mag * nij[1],
            v_n_mag * nij[2]
        };

        double v_t[3] = {
            vij[0] - v_n[0],
            vij[1] - v_n[1],
            vij[2] - v_n[2]
        };

        // --------------------------------
        //  lubrication force 
        // --------------------------------
        double F_n[3] = {0.0, 0.0, 0.0};

        // Translational resistance
        F_n[0] += (X11A * v_n[0] + Y11A * v_t[0]) * -1;
        F_n[1] += (X11A * v_n[1] + Y11A * v_t[1]) * -1;
        F_n[2] += (X11A * v_n[2] + Y11A * v_t[2]) * -1;

        // wi × nij
        double wixn[3] = {
            wi[1]*nij[2] - wi[2]*nij[1],
            wi[2]*nij[0] - wi[0]*nij[2],
            wi[0]*nij[1] - wi[1]*nij[0]
        };

        // wj × nij
        double wjxn[3] = {
            wj[1]*nij[2] - wj[2]*nij[1],
            wj[2]*nij[0] - wj[0]*nij[2],
            wj[0]*nij[1] - wj[1]*nij[0]
        };

        F_n[0] += Y11B * wixn[0] + Y21B * wjxn[0];
        F_n[1] += Y11B * wixn[1] + Y21B * wjxn[1];
        F_n[2] += Y11B * wixn[2] + Y21B * wjxn[2];

        // Multiply by viscosity
        F_n[0] *= eta;
        F_n[1] *= eta;
        F_n[2] *= eta;

        // --------------------------------
        // lubrication torque 
        // --------------------------------
        double F_t[3] = {0.0, 0.0, 0.0};

        // vij × nij
        double vixn[3] = {
            vij[1]*nij[2] - vij[2]*nij[1],
            vij[2]*nij[0] - vij[0]*nij[2],
            vij[0]*nij[1] - vij[1]*nij[0]
        };

        F_t[0] += Y11B * vixn[0] * -1;
        F_t[1] += Y11B * vixn[1]* -1;
        F_t[2] += Y11B * vixn[2]* -1;

        // Rotational term: (I - nn)(Y11C wi + Y12C wj)
        double omega_term[3] = {
            Y11C * wi[0] + Y12C * wj[0],
            Y11C * wi[1] + Y12C * wj[1],
            Y11C * wi[2] + Y12C * wj[2]
        };

        double omega_n =
            omega_term[0]*nij[0] +
            omega_term[1]*nij[1] +
            omega_term[2]*nij[2];

        double omega_t[3] = {
            omega_term[0] - omega_n * nij[0],
            omega_term[1] - omega_n * nij[1],
            omega_term[2] - omega_n * nij[2]
        };

        F_t[0] -= omega_t[0];
        F_t[1] -= omega_t[1];
        F_t[2] -= omega_t[2];

        // Multiply by viscosity
        F_t[0] *= eta;
        F_t[1] *= eta;
        F_t[2] *= eta;

        // --------------------------------
        // Total lubrication force
        // --------------------------------
        F_lub[0] = F_n[0];
        F_lub[1] = F_n[1];
        F_lub[2] = F_n[2];
        t_lub[0] = F_t[0];
        t_lub[1] = F_t[1];
        t_lub[2] = F_t[2];
    }

    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {
      const int i = scdata.i;
      const int j = scdata.j;
      const int itype = scdata.itype;
      const int jtype = scdata.jtype;
      const double radi = scdata.radi;
      const double radj = scdata.is_wall ? radi : scdata.radj;
      const double r = sqrt(scdata.rsq);
      const double radsum = scdata.radsum;
      const double dist = scdata.is_wall ? r - radi : r - (radi + radj);
      const double rEff = radi*radj / (radi+radj);

        double const *x_i = atom->x[i];
        double const *x_j = atom->x[j];

        double const *omega_i = atom->omega[i];
        double const *omega_j = atom->omega[j];

        const double dx = scdata.delta[0];
        const double dy = scdata.delta[1];
        const double dz = scdata.delta[2];
         const double *v_i = scdata.v_i;
        const double  *v_j = scdata.v_j;

        double eta=fluidDynamicViscosity; //WATER
        double h = sqrt(dx*dx+dy*dy+dz*dz) -(radi+radj);
        double epsi= 2*h/(radj+radi);
        double beta=radj/radi;
        double F_lub[3];
        double T_lub[3];
        computeLubricationForce(radi,radj,x_i,x_j,v_i,v_j,omega_i,omega_j,eta,epsi,beta,F_lub,T_lub);
        scdata.has_force_update = true;
        if(scdata.is_wall) {
            const double area_ratio = scdata.area_ratio;
            i_forces.delta_F[0] += F_lub[0] * area_ratio;
            i_forces.delta_F[1] += F_lub[1] * area_ratio;
            i_forces.delta_F[2] += F_lub[2] * area_ratio;
            i_forces.delta_torque[0] += T_lub[0] * area_ratio;
            i_forces.delta_torque[1] += T_lub[1] * area_ratio;
            i_forces.delta_torque[2] += T_lub[2] * area_ratio;
        } else {
            i_forces.delta_F[0] += F_lub[0];
            i_forces.delta_F[1] += F_lub[1];
            i_forces.delta_F[2] += F_lub[2];
            i_forces.delta_torque[0] += T_lub[0]; 
            i_forces.delta_torque[1] += T_lub[1];
            i_forces.delta_torque[2] += T_lub[2];

            j_forces.delta_F[0] -= F_lub[0];
            j_forces.delta_F[1] -= F_lub[1];
            j_forces.delta_F[2] -= F_lub[2];
            j_forces.delta_torque[0] -= T_lub[0]; 
            j_forces.delta_torque[1] -= T_lub[1];
            j_forces.delta_torque[2] -= T_lub[2];
        }
    }
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double ** k_n;
    double ** k_t;
    double ** gamma_n;
    double ** gamma_t;
    double fluidDynamicViscosity;

    bool tangential_damping;
    bool limitForce;
    bool displayedSettings;
  };
}
}
#endif // NORMAL_MODEL_HERTZ_STIFFNESS_H_
#endif
