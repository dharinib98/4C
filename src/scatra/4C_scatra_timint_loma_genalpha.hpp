// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_TIMINT_LOMA_GENALPHA_HPP
#define FOUR_C_SCATRA_TIMINT_LOMA_GENALPHA_HPP

#include "4C_config.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_scatra_timint_genalpha.hpp"
#include "4C_scatra_timint_loma.hpp"

FOUR_C_NAMESPACE_OPEN


namespace ScaTra
{
  class TimIntLomaGenAlpha : public virtual ScaTraTimIntLoma, public virtual TimIntGenAlpha
  {
   public:
    /// Standard Constructor
    TimIntLomaGenAlpha(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);


    /// initialize time integration scheme
    void init() override;

    /// setup time integration scheme
    void setup() override;

    /// predict thermodynamic pressure and time derivative for low-Mach-number flow
    void predict_therm_pressure() override;

    /// compute values of thermodynamic pressure at intermediate time steps
    void compute_therm_pressure_intermediate_values() override;

    /// compute thermodynamic pressure and time derivative for low-Mach-number flow
    void compute_therm_pressure() override;

    ///  compute time derivative of thermodynamic pressure
    void compute_therm_pressure_time_derivative() override;

    /// update thermodynamic pressure and time derivative for low-Mach-number flow
    void update_therm_pressure() override;

    /// read restart data
    void read_restart(
        const int step, std::shared_ptr<Core::IO::InputControl> input = nullptr) override;

    /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
    double therm_press_af() override { return thermpressaf_; }

    /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
    double therm_press_am() override { return thermpressam_; }

    /// routine to return time derivative of thermo. press. at time step n+alpha_F for
    /// low-Mach-number flow
    double therm_press_dt_af() override { return thermpressdtaf_; }

    /// routine to return time derivative of thermo. press. at time step n+alpha_M for
    /// low-Mach-number flow
    double therm_press_dt_am() override { return thermpressdtam_; }

   protected:
    void write_restart() const override;

    /// dynamic Smagorinsky model
    void dynamic_computation_of_cs() override;

    /// dynamic Vreman model
    void dynamic_computation_of_cv() override;

    /// add thermodynamic pressure to parameter list for element evaluation
    void add_therm_press_to_parameter_list(Teuchos::ParameterList& params  //!< parameter list
        ) override;

   private:
    /// LOMA-specific parameter: thermodynamic pressure at n+alpha_F and n+alpha_M
    /// and its time derivative at n+alpha_F and n+alpha_M
    double thermpressaf_;
    double thermpressam_;
    double thermpressdtaf_;
    double thermpressdtam_;


  };  // class TimIntLomaGenAlpha

}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
