/*
 *  stdp_tsodyks_connection.h
 *
 *  This is a customized synapse module merging stdp_synapse
 *  and tsodyks_synapse.
 *
 *  This custom module is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef STDP_TSODYKS_CONNECTION_H
#define STDP_TSODYKS_CONNECTION_H

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"

namespace nest
{

/** @BeginDocumentation
@ingroup Synapses
@ingroup stdp

Name: stdp_tsodyks_synapse - Synapse type for spike-timing dependent
plasticity and Tsodyks short-term plasticity.

Description:

This module is a simple combinatin of stdp_connection.h and 
tsodyks_connection.h. For detail, please check those individual 
files. The long-term plasticity example are:

Examples:

    multiplicative STDP [2]  mu_plus = mu_minus = 1.0
    additive STDP       [3]  mu_plus = mu_minus = 0.0
    Guetig STDP         [1]  mu_plus = mu_minus = [0.0,1.0]
    van Rossum STDP     [4]  mu_plus = 0.0 mu_minus = 1.0

Parameters:
\verbatim embed:rst
========= =======  ======================================================
 tau_plus  ms      Time constant of STDP window, potentiation
                   (tau_minus defined in post-synaptic neuron)
 lambda    real    Step size
 alpha     real    Asymmetry parameter (scales depressing increments as
                   alpha*lambda)
 mu_plus   real    Weight dependence exponent, potentiation
 mu_minus  real    Weight dependence exponent, depression
 Wmax      real    Maximum allowed weight
 U         real    Maximum probability of release [0,1]
 tau_psc   ms      Time constant of synaptic current
 tau_fac   ms      Time constant for facilitation
 tau_rec   ms      Time constant for depression
 x         real    Initial fraction of synaptic vesicles in the readily
                   releasable pool [0,1]
 y         real    Initial fraction of synaptic vesicles in the synaptic
                   cleft [0,1]
========= =======  ======================================================
\endverbatim

Transmits: SpikeEvent

References:

\verbatim embed:rst
.. [1] Guetig et al. (2003). Learning input correlations through nonlinear
       temporally asymmetric hebbian plasticity. Journal of Neuroscience,
       23:3697-3714 DOI: https://doi.org/10.1523/JNEUROSCI.23-09-03697.2003
.. [2] Rubin J, Lee D, Sompolinsky H (2001). Equilibrium
       properties of temporally asymmetric Hebbian plasticity. Physical Review
       Letters, 86:364-367. DOI: https://doi.org/10.1103/PhysRevLett.86.364
.. [3] Song S, Miller KD, Abbott LF (2000). Competitive Hebbian learning
       through spike-timing-dependent synaptic plasticity. Nature Neuroscience
       3(9):919-926.
       DOI: https://doi.org/10.1038/78829
.. [4] van Rossum MCW, Bi G-Q, Turrigiano GG (2000). Stable Hebbian learning
       from spike timing-dependent plasticity. Journal of Neuroscience,
       20(23):8812-8821.
       DOI: https://doi.org/10.1523/JNEUROSCI.20-23-08812.2000
\endverbatim

FirstVersion: March 2006

Author: Moritz Helias, Abigail Morrison

Adapted by: Philipp Weidel

SeeAlso: synapsedict, tsodyks_synapse, static_synapse
*/

// connections are templates of target identifier type (used for pointer /
// target index addressing) derived from generic connection template
template < typename targetidentifierT >
class STDPTsodyksConnection : public Connection< targetidentifierT >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDPTsodyksConnection();

  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDPTsodyksConnection( const STDPTsodyksConnection& );

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param cp common properties of all synapses (empty).
   */
  void send( Event& e, thread t, const CommonSynapseProperties& cp );

  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport )
    {
      return invalid_port_;
    }
  };

  void
  check_connection( Node& s, Node& t, rport receptor_type, const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;

    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike_ - get_delay(), get_delay() );
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }

private:
  double
  facilitate_( double w, double kplus )
  {
    double norm_w = ( w / Wmax_ ) + ( lambda_ * std::pow( 1.0 - ( w / Wmax_ ), mu_plus_ ) * kplus );
    return norm_w < 1.0 ? norm_w * Wmax_ : Wmax_;
  }

  double
  depress_( double w, double kminus )
  {
    double norm_w = ( w / Wmax_ ) - ( alpha_ * lambda_ * std::pow( w / Wmax_, mu_minus_ ) * kminus );
    return norm_w > 0.0 ? norm_w * Wmax_ : 0.0;
  }

  // data members of each connection
  double weight_;
  double tau_plus_;
  double lambda_;
  double alpha_;
  double mu_plus_;
  double mu_minus_;
  double Wmax_;
  double Kplus_;

  double t_lastspike_;

  /*
  The following lines of code are copied from tsodyks_synapse.h 
  */
  
  double tau_psc_;     //!< [ms] time constant of postsyn current
  double tau_fac_;     //!< [ms] time constant for fascilitation
  double tau_rec_;     //!< [ms] time constant for recovery
  double U_;           //!< asymptotic value of probability of release
  double x_;           //!< amount of resources in recovered state
  double y_;           //!< amount of resources in active state
  double u_;           //!< actual probability of release
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param t The thread on which this connection is stored.
 * \param cp Common properties object, containing the stdp parameters.
 */
template < typename targetidentifierT >
inline void
STDPTsodyksConnection< targetidentifierT >::send( Event& e, thread t, const CommonSynapseProperties& )
{
  // synapse STDP depressing/facilitation dynamics
  const double t_spike = e.get_stamp().get_ms();
  const double h = t_spike - t_lastspike_;  // From tsodyks_synapse.h

  // use accessor functions (inherited from Connection< >) to obtain delay and
  // target
  Node* target = get_target( t );
  double dendritic_delay = get_delay();

  // !!! Tsodyks Model

  // t_lastspike_ = 0 initially
  // this has no influence on the dynamics, IF y = z = 0 initially
  // !!! x != 1.0 -> z != 0.0 -> t_lastspike_=0 has influence on dynamics

  /*
  The following lines of code are dynamics of the Tsodyks's model 
  copied from tsodyks_synapse.h 
  */

  // propagator
  // TODO: use expm1 here instead, where applicable
  double Puu = ( tau_fac_ == 0.0 ) ? 0.0 : std::exp( -h / tau_fac_ );
  double Pyy = std::exp( -h / tau_psc_ );
  double Pzz = std::exp( -h / tau_rec_ );

  double Pxy = ( ( Pzz - 1.0 ) * tau_rec_ - ( Pyy - 1.0 ) * tau_psc_ ) / ( tau_psc_ - tau_rec_ );
  double Pxz = 1.0 - Pzz;

  double z = 1.0 - x_ - y_;

  // propagation t_lastspike_ -> t_spike
  // don't change the order !

  u_ *= Puu;
  x_ += Pxy * y_ + Pxz * z;
  y_ *= Pyy;

  // delta function u
  u_ += U_ * ( 1.0 - u_ );

  // postsynaptic current step caused by incoming spike
  double delta_y_tsp = u_ * x_;

  // delta function x, y
  x_ -= delta_y_tsp;
  y_ += delta_y_tsp;

  /*
	The following is from stdp_connection.h
  */

  // get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque< histentry >::iterator start;
  std::deque< histentry >::iterator finish;

  double old_weight_ = weight_;		// This variable is used to store the pre-stdp weight

  // For a new synapse, t_lastspike_ contains the point in time of the last
  // spike. So we initially read the
  // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of
  // history[0, ..., t_last_spike - dendritic_delay] have been
  // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for
  // details.
  target->get_history( t_lastspike_ - dendritic_delay, t_spike - dendritic_delay, &start, &finish );
  // facilitation due to post-synaptic spikes since last pre-synaptic spike
  double minus_dt;
  while ( start != finish )
  {
    minus_dt = t_lastspike_ - ( start->t_ + dendritic_delay );
    ++start;
    // get_history() should make sure that
    // start->t_ > t_lastspike - dendritic_delay, i.e. minus_dt < 0
    assert( minus_dt < -1.0 * kernel().connection_manager.get_stdp_eps() );
    weight_ = facilitate_( weight_, Kplus_ * std::exp( minus_dt / tau_plus_ ) );
  }

  const double _K_value = target->get_K_value( t_spike - dendritic_delay );
  weight_ = depress_( weight_, _K_value );

  e.set_receiver( *target );
  e.set_weight( delta_y_tsp * old_weight_ ); // Changed for short-term plasticity
  // use accessor functions (inherited from Connection< >) to obtain delay in
  // steps and rport
  e.set_delay_steps( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

  Kplus_ = Kplus_ * std::exp( ( t_lastspike_ - t_spike ) / tau_plus_ ) + 1.0;

  t_lastspike_ = t_spike;
}


template < typename targetidentifierT >
STDPTsodyksConnection< targetidentifierT >::STDPTsodyksConnection()
  : ConnectionBase()
  , weight_( 1.0 )
  , tau_plus_( 20.0 )
  , lambda_( 0.01 )
  , alpha_( 1.0 )
  , mu_plus_( 1.0 )
  , mu_minus_( 1.0 )
  , Wmax_( 100.0 )
  , Kplus_( 0.0 )
  // Head of the block from tsodyks_connection.h
  , tau_psc_( 3.0 )
  , tau_fac_( 0.0 )
  , tau_rec_( 800.0 )
  , U_( 0.5 )
  , x_( 1.0 )
  , y_( 0.0 )
  , u_( 0.0 )
  // End of the block from tsodyks_connection.h
  , t_lastspike_( 0.0 )
{
}

template < typename targetidentifierT >
STDPTsodyksConnection< targetidentifierT >::STDPTsodyksConnection( const STDPTsodyksConnection< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , tau_plus_( rhs.tau_plus_ )
  , lambda_( rhs.lambda_ )
  , alpha_( rhs.alpha_ )
  , mu_plus_( rhs.mu_plus_ )
  , mu_minus_( rhs.mu_minus_ )
  , Wmax_( rhs.Wmax_ )
  , Kplus_( rhs.Kplus_ )
  // Head of the block from tsodyks_connection.h
  , tau_psc_( rhs.tau_psc_ )
  , tau_fac_( rhs.tau_fac_ )
  , tau_rec_( rhs.tau_rec_ )
  , U_( rhs.U_ )
  , x_( rhs.x_ )
  , y_( rhs.y_ )
  , u_( rhs.u_ )
  // End of the block from tsodyks_connection.h
  , t_lastspike_( rhs.t_lastspike_ )
{
}

template < typename targetidentifierT >
void
STDPTsodyksConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, names::tau_plus, tau_plus_ );
  def< double >( d, names::lambda, lambda_ );
  def< double >( d, names::alpha, alpha_ );
  def< double >( d, names::mu_plus, mu_plus_ );
  def< double >( d, names::mu_minus, mu_minus_ );
  def< double >( d, names::Wmax, Wmax_ );

  // Head of the block from tsodyks_connection.h
  def< double >( d, names::U, U_ );
  def< double >( d, names::tau_psc, tau_psc_ );
  def< double >( d, names::tau_rec, tau_rec_ );
  def< double >( d, names::tau_fac, tau_fac_ );
  def< double >( d, names::x, x_ );
  def< double >( d, names::y, y_ );
  def< double >( d, names::u, u_ );
  // End of the block from tsodyks_connection.h
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
STDPTsodyksConnection< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  // Head of the block from tsodyks_connection.h

  // Handle parameters that may throw an exception first, so we can leave the
  // synapse untouched
  // in case of invalid parameter values
  double x = x_;
  double y = y_;
  updateValue< double >( d, names::x, x );
  updateValue< double >( d, names::y, y );

  if ( x + y > 1.0 )
  {
    throw BadProperty( "x + y must be <= 1.0." );
  }

  x_ = x;
  y_ = y;

  // End of the block from tsodyks_connection.h

  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::tau_plus, tau_plus_ );
  updateValue< double >( d, names::lambda, lambda_ );
  updateValue< double >( d, names::alpha, alpha_ );
  updateValue< double >( d, names::mu_plus, mu_plus_ );
  updateValue< double >( d, names::mu_minus, mu_minus_ );
  updateValue< double >( d, names::Wmax, Wmax_ );

  // Head of the block from tsodyks_connection.h
  updateValue< double >( d, names::U, U_ );
  if ( U_ > 1.0 || U_ < 0.0 )
  {
    throw BadProperty( "U must be in [0,1]." );
  }

  updateValue< double >( d, names::tau_psc, tau_psc_ );
  if ( tau_psc_ <= 0.0 )
  {
    throw BadProperty( "tau_psc must be > 0." );
  }

  updateValue< double >( d, names::tau_rec, tau_rec_ );
  if ( tau_rec_ <= 0.0 )
  {
    throw BadProperty( "tau_rec must be > 0." );
  }

  updateValue< double >( d, names::tau_fac, tau_fac_ );
  if ( tau_fac_ < 0.0 )
  {
    throw BadProperty( "tau_fac must be >= 0." );
  }

  updateValue< double >( d, names::u, u_ );

  // End of the block from tsodyks_connection.h

  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) ) == ( ( Wmax_ >= 0 ) - ( Wmax_ < 0 ) ) ) )
  {
    throw BadProperty( "Weight and Wmax must have same sign." );
  }
}

} // of namespace nest

#endif // of #ifndef STDP_TSODYKS_CONNECTION_H
