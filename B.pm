package Physics::UEMColumn::ThreeQuartersButlerAccelerator;

use Moose;
use namespace::autoclean;

BEGIN { extends 'Physics::UEMColumn::Accelerator'; }

use Method::Signatures;

use Physics::UEMColumn::Auxiliary ':all';
use Math::Trig qw/tanh sech/;
use MooseX::Types::NumUnit qw/num_of_unit/;

has 'voltage' => ( isa => num_of_unit('V'), is => 'ro', required => 1 );
has '+length' => ( isa => num_of_unit('m'), is => 'ro', required => 0, default =>  0.04570984 . 'm' );
has 'sharpness' => ( isa => 'Num', is => 'ro', default => 1 );


# Cyl gun fit
my $pm_a = 0.544667;
my $pm_b = 0.00278104;
my $pm_c = 332.551;
my $pm_d = 0.0415734;
my $pm_e = 213.238;
my $pm_f = 0.0168457;
my $pm_g = 1761.42;
my $pm_h = 31.4205;

if (1) {
   $pm_a = 0.7323682193121402;
   $pm_b = 0.0018258141203642293;
   $pm_c = 419.22305621512435;
   $pm_d = 0.042697578590972284;
   $pm_e = 231.22817596631087;
   $pm_f = 0.016856140662231493;
   $pm_g = 1259.5318887451284;
   $pm_h = 28.65;
}

method field () {
  my $acc_length = 1/$pm_h;
  return $self->voltage/$acc_length;

  my $pulse_z = 0;
  return $self->voltage * $pm_h * ( ($pm_a*(1 - tanh($pm_c*(-$pm_b + $pulse_z))))/2 +  ((1 + tanh($pm_c*(-$pm_b + $pulse_z)))*(1 - tanh($pm_e*(-$pm_d + $pulse_z))))/  (4*exp($pm_g*(-$pm_f + $pulse_z)**2)) )
}

method effect () {
  my $anode_pos = $self->length;
  my $acc_voltage = $self->voltage;
  my $force = qe * $acc_voltage;

  # cutoff is used oddly here
  my $cutoff = $self->cutoff;

  my $acc = sub {
    my ($t, $pulse_z, $pulse_v) = @_;
    if ($pulse_z / $anode_pos > $cutoff) {

      return 0;
    }
    return $force / me * $pm_h * ( ($pm_a*(1 - tanh($pm_c*(-$pm_b + $pulse_z))))/2 +  ((1 + tanh($pm_c*(-$pm_b + $pulse_z)))*(1 - tanh($pm_e*(-$pm_d + $pulse_z))))/  (4*exp($pm_g*(-$pm_f + $pulse_z)**2)) );

  };

  my $acc_mt = sub {
    my ($t, $pulse_z, $pulse_v) = @_;

    if ($pulse_z / $anode_pos > $cutoff) {
      return 0;
    }

    return $force/2 * $pm_h * ((-1 - tanh($pm_c*(-$pm_b +$pulse_z)))*($pm_e*sech($pm_e*(-$pm_d +$pulse_z))**2 + 2*$pm_g*($pm_f -$pulse_z)*(-1 + tanh($pm_e*(-$pm_d +$pulse_z)))) - $pm_c*sech($pm_c*(-$pm_b +$pulse_z))**2*(-1 + 2*exp($pm_g*($pm_f -$pulse_z)**2)*$pm_a + tanh($pm_e*(-$pm_d +$pulse_z))))/(4*exp($pm_g*($pm_f -$pulse_z)**2));
  };

  my $acc_mz = sub {
    my ($t, $pulse_z, $pulse_v) = @_;

    if ($pulse_z / $anode_pos > $cutoff) {
      return 0;
    }

    return -$force * $pm_h * ((-1 - tanh($pm_c*(-$pm_b +$pulse_z)))*($pm_e*sech($pm_e*(-$pm_d +$pulse_z))**2 + 2*$pm_g*($pm_f -$pulse_z)*(-1 + tanh($pm_e*(-$pm_d +$pulse_z)))) - $pm_c*sech($pm_c*(-$pm_b +$pulse_z))**2*(-1 + 2*exp($pm_g*($pm_f -$pulse_z)**2)*$pm_a + tanh($pm_e*(-$pm_d +$pulse_z))))/(4*exp($pm_g*($pm_f -$pulse_z)**2));
  };

  #TODO add anode effects
  return {acc => $acc, M_t => $acc_mt, M_z => $acc_mz};

}

method est_exit_vel () {
  return sqrt( 2 * qe * $self->voltage / me );
}

method est_exit_time () {
  # assumes pulse has initial vel zero
  return $self->length() * sqrt( 2 * me / ( qe * $self->voltage ) );
}

__PACKAGE__->meta->make_immutable;

1;

=head1 NAME

Physics::UEMColumn::DCAccelerator - A class representing a DC acceleration region in a UEM system

=head1 SYNOPSIS

 use Physics::UEMColumn alias => ':standard';
 my $acc = DCAccelerator->new(
   length  => '20 mm',
   voltage => '20 kilovolts',
 );

=head1 DESCRIPTION

L<Physics::UEMColumn::Accelerator> is a class representing a DC (static electric field) acceleration region in a UEM system. It is a subclass of L<Physics::UEMColumn::Accelerator> and inherits its attributes and methods. Additionally it provides:

=head1 ATTRIBUTES

=over

=item C<voltage>

The static electric potential in the accelerator. Unit: V

=item C<sharpness>

The potential is modeled as a C<tanh>, this parameter (defaults to 10) is related to the slope of the tanh near the end of the region. For example a value approaching infinity approximates a step function.

=back

=head1 METHODS

=over

=item C<field>

Defined as C<voltage> / C<length>

=item C<effect>

Returns a hash reference of effect subroutine references (C<M_t>, C<M_z>, C<acc_z>). See L<Physics::UEMColumn::Element/METHODS> for more.

=item C<est_exit_vel>

Returns an estimate of the velocity of the pulse on exiting the region. This in not likely to be exact. It is used in estimating the end time of the simulation. This overrides the base class and is specific to DC accelerators.

=item C<est_exit_time>

Returns an estimate of the time that the pulse on exits the region. This in not likely to be exact. It is used in estimating the end time of the simulation. This overrides the base class and is specific to DC accelerators.

=back

=head1 SOURCE REPOSITORY

L<http://github.com/jberger/Physics-UEMColumn>

=head1 AUTHOR

Joel Berger, E<lt>joel.a.berger@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012-2013 by Joel Berger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

