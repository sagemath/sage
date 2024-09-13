"""
Units of measurement

This is the units package. It contains information about many units
and conversions between them.

TUTORIAL:

To return a unit::

    sage: units.length.meter
    meter

This unit acts exactly like a symbolic variable::

    sage: s = units.length.meter
    sage: s^2
    meter^2
    sage: s + var('x')
    meter + x

Units have additional information in their docstring::

    sage: # You would type: units.force.dyne?
    sage: print(units.force.dyne.__doc__)
    CGS unit for force defined to be gram*centimeter/second^2.
    Equal to 10^-5 newtons.


You may call the convert function with units::

    sage: t = units.mass.gram*units.length.centimeter/units.time.second^2
    sage: t.convert(units.mass.pound*units.length.foot/units.time.hour^2)
    5400000000000/5760623099*(foot*pound/hour^2)
    sage: t.convert(units.force.newton)
    1/100000*newton

Calling the convert function with no target returns base SI units::

    sage: t.convert()
    1/100000*kilogram*meter/second^2

Giving improper units to convert to raises a :exc:`ValueError`::

    sage: t.convert(units.charge.coulomb)
    Traceback (most recent call last):
    ...
    ValueError: Incompatible units

Converting temperatures works as well::

    sage: s = 68*units.temperature.fahrenheit
    sage: s.convert(units.temperature.celsius)
    20*celsius
    sage: s.convert()
    293.150000000000*kelvin

Trying to multiply temperatures by another unit then converting raises a :exc:`ValueError`::

    sage: wrong = 50*units.temperature.celsius*units.length.foot
    sage: wrong.convert()
    Traceback (most recent call last):
    ...
    ValueError: cannot convert

TESTS:

Check that :issue:`12373` is fixed::

    sage: b = units.amount_of_substance.mole
    sage: b.convert(units.amount_of_substance.elementary_entity)
    6.02214129000000e23*elementary_entity

AUTHORS:

- David Ackerman
- William Stein
"""

###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2009 David Ackerman <davidnackerman@gmail.com>
#                          William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
###############################################################################

# standard Python libraries
import re

# Sage library
from .ring import SR
from .expression import Expression
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.misc.instancedoc import instancedoc
from sage.rings.rational_field import QQ

###############################################################################
# Unit conversions dictionary.
###############################################################################

unitdict = {
'acceleration':
        {'gal': '1/100',
        'galileo': '1/100',
        'gravity': '9.80665000000000'},

'amount_of_substance':
        {'elementary_entity': '1/6.02214129000000e23',
        'mole': 1},

'angles':
        {'arc_minute': '1/10800*pi',
        'arc_second': '1/648000*pi',
        'degree': '1/180*pi',
        'grade': '1/200*pi',
        'quadrant': '1/2*pi',
        'radian': 1,
        'right_angle': '1/2*pi'},

'area':
        {'acre': '316160658/78125',
        'are': 100,
        'barn': '1/10000000000000000000000000000',
        'hectare': 10000,
        'rood': '158080329/156250',
        'section': '40468564224/15625',
        'square_chain': '158080329/390625',
        'square_meter': 1,
        'township': '1456868312064/15625'},

'capacitance':
        {'abfarad': 1000000000,
        'farad': 1,
        'statfarad': '25000/22468879468420441'},

'charge':
        {'abcoulomb': 10,
        'coulomb': 1,
        'elementary_charge': '1.60217646200000e-19',
        'faraday': '96485.3399000000',
        'franklin': '1/2997924580',
        'statcoulomb': '1/2997924580'},

'conductance':
        {'abmho': 1000000000,
        'mho': 1,
        'siemens': 1},

'current':
        {'abampere': 10,
        'amp': 1,
        'ampere': 1,
        'biot': 10,
        'statampere': '1/2997924580'},

'electric_potential':
        {'abvolt': '1/100000000',
        'statvolt': '149896229/500000',
        'volt': 1},

'energy':
        {'british_thermal_unit': '52752792631/50000000',
        'btu': '52752792631/50000000',
        'calorie': '10467/2500',
        'electron_volt': '1.60217733000000e-19',
        'erg': '1/10000000',
        'ev': '1.60217733000000e-19',
        'joule': 1,
        'rydberg': '2.17987200000000e-18',
        'therm': '52752792631/500'},

'fiber_linear_mass_density':
        {'denier': '1/9000000',
        'tex': '1/1000000'},

'force':
        {'dyne': '1/100000',
        'gram_weight': '196133/20000000',
        'kilogram_force': '196133/20000',
        'kilogram_weight': '196133/20000',
        'newton': 1,
        'pound_force': '8896443230521/2000000000000',
        'pound_weight': '8896443230521/2000000000000',
        'poundal': '17281869297/125000000000',
        'ton_force': '8896443230521/1000000000'},

'frequency':
        {'1/second': 1,
        'hertz': 1},

'illuminance':
        {'foot_candle': '1562500/145161',
        'lux': 1,
        'phot': 10000},

'inductance':
        {'abhenry': '1/1000000000',
        'henry': 1,
        'stathenry': '22468879468420441/25000'},

'information':
        {'bit': 1,
        'byte': 8,
        'nibble': 4},

'information_rate':
        {'baud': 1},

'inverse_length':
        {'diopter': 1,
        'kayser': 100},

'length':
        {'angstrom': '1/10000000000',
        'astronomical_unit': '149597870691',
        'bolt': '4572/125',
        'cable_international': '926/5',
        'cable_us': '27432/125',
        'caliber': '127/500000',
        'centimeter': '1/100',
        'chain': '12573/625',
        'cicero': '125/27706',
        'cubit': '1143/2500',
        'didot': '125/332472',
        'dtp_point': '127/360000',
        'ell': '1143/1000',
        'fathom': '1143/625',
        'feet': '381/1250',
        'fermi': '1/1000000000000000',
        'foot': '381/1250',
        'furlong': '25146/125',
        'hand': '127/1250',
        'inch': '127/5000',
        'kilometer': 1000,
        'league': '603504/125',
        'light_year': '9460730472580800',
        'link': '12573/62500',
        'meter': 1,
        'micron': '1/1000000',
        'mil': '127/5000000',
        'millimeter': '1/1000',
        'mile': '201168/125',
        'nautical_mile': 1852,
        'parsec': '3.08570000000000e16',
        'perch': '12573/2500',
        'pica': '127/30000',
        'pole': '12573/2500',
        'rod': '12573/2500',
        'rope': '762/125',
        'skein': '13716/125',
        'stadion': '118491/625',
        'stadium': '115443/625',
        'statute_mile': '201168/125',
        'survey_foot': '1200/3937',
        'survey_mile': '6336000/3937',
        'x_unit': '1.00210000000000e-13',
        'yard': '1143/1250'},

'luminance':
        {'apostilb': '1/pi',
        'lambert': '10000/pi',
        'nit': 1,
        'stilb': 10000},

'luminous_energy':
        {'lumerg': 1,
        'talbot': 1},

'luminous_flux':
        {'lumen': 1},

'luminous_intensity':
        {'candela': 1,
        'candle': 1,
        'hefnerkerze': '1019/1128'},

'magnetic_field':
        {'gauss': '1/10000',
        'tesla': 1},

'magnetic_flux':
        {'maxwell': '1/100000000',
        'weber': 1},

'magnetic_intensity':
        {'oersted': '250/pi'},

'magnetic_moment':
        {'bohr_magneton': '9.27400915000000e-24',
        'nuclear_magneton': '5.05078324000000e-27'},

'magnetomotive_force':
        {'ampere_turn': 1,
        'gilbert': '5/2/pi'},

'mass':
        {'amu': '1.66053878200000e-27',
        'assay_ton': '7/240',
        'atomic_mass_unit': '1.66053878200000e-27',
        'avoirdupois_ounce': '45359237/1600000000',
        'avoirdupois_pound': '45359237/100000000',
        'bale': '45359237/200000',
        'carat': '1/5000',
        'cental': '45359237/1000000',
        'dalton': '1.66053878200000e-27',
        'drachma': "(0.00429234000000000, {'greek':1})",
        'geepound': '14593903/1000000',
        'grain': '6479891/100000000000',
        'gram': '1/1000',
        'gross_hundredweight': '317514659/6250000',
        'hundredweight': '317514659/6250000',
        'kilogram': 1,
        'libra': '0.325971000000000',
        'long_ton': '317514659/312500',
        'metric_ton': '1000',
        'mina': "(0.429234000000000, {'greek':100})",
        'net_hundredweight': '45359237/1000000',
        'obol': "(0.000715380000000000,{'greek':1/6})",
        'ounce': '45359237/1600000000',
        'ounce_troy': '19439673/625000000',
        'pennyweight': '19439673/12500000000',
        'pondus': '0.325969000000000',
        'pound': '45359237/100000000',
        'pound_troy': '58319019/156250000',
        'quintal': '100',
        'shekel': '0.0141000000000000',
        'short_hundredweight': '45359237/1000000',
        'short_ton': '45359237/50000',
        'slug': '14593903/1000000',
        'solar_mass': '1.98892000000000e30',
        'stone': '317514659/50000000',
        'talent': "(25.7540400000000, {'greek':6000})",
        'ton': '45359237/50000',
        'tonne': '1000',
        'wey': '2857631931/25000000'},

'power':
        {'cheval_vapeur': '588399/800',
        'horsepower': '37284993579113511/50000000000000',
        'watt': 1},

'pressure':
        {'atmosphere': 101325,
        'bar': 100000,
        'barye': QQ((1, 10)),
        'inch_mercury': '3386.38900000000',
        'millimeter_mercury': '133.322400000000',
        'mmhg': '133.322400000000',
        'pa': 1,
        'pascal': 1,
        'pounds_per_square_inch': '8896443230521/1290320000',
        'psi': '8896443230521/1290320000',
        'torr': '20265/152'},

'radiation':
        {'becquerel': 1,
        'curie': 37000000000,
        'rutherford': 1000000},

'radiation_absorbed':
        {'gray': 1,
        'rad': '1/100'},

'radiation_ionizing':
        {'roentgen': '0.000258000000000000',
        'rontgen': '0.000258000000000000'},

'resistance':
        {'abohm': '1/1000000000',
        'ohm': 1,
        'statohm': '22468879468420441/25000'},

'si_prefixes':
        {'atto': '1/1000000000000000000',
        'centi': '1/100',
        'deca': 10,
        'deci': QQ((1, 10)),
        'exa': '1000000000000000000',
        'femto': '1/1000000000000000',
        'giga': '1000000000',
        'hecto': 100,
        'kilo': 1000,
        'mega': 1000000,
        'micro': '1/1000000',
        'milli': '1/1000',
        'nano': '1/1000000000',
        'peta': '1000000000000000',
        'pico': '1/1000000000000',
        'tera': 1000000000000,
        'yocto': '1/1000000000000000000000000',
        'yotta': '1000000000000000000000000',
        'zepto': '1/1000000000000000000000',
        'zetta': '1000000000000000000000'},

'solid_angle':
        {'steradian': 1},

'temperature':
        {'celsius': '(x + 273.15), (x), (x*9/5 + 32), ((x+273.15)*9/5)',
        'centigrade': '(x + 273.15), (x), (x*9/5 + 32), ((x+273.15)*9/5)',
        'fahrenheit': '(5/9*(x + 459.67)), ((x - 32)*5/9), (x), (x+459.67)',
        'kelvin': '(x), (x - 273.15), (x*9/5 - 459.67), (x*9/5)',
        'rankine': '(5/9*x), ((x-491.67)*5/9), (x-459.67), (x)'},

'time':
        {'century': 3153600000,
        'day': 86400,
        'decade': 315360000,
        'fortnight': 1209600,
        'hour': 3600,
        'millenium': 31536000000,
        'minute': 60,
        'month': 2628000,
        'second': 1,
        'sidereal_day': "(86164.0905308330, {'sidereal':86400})",
        'sidereal_second': "(0.997269566329086, {'sidereal':1})",
        'sidereal_year': '3.15581497632000e7',
        'tropical_year': '3.15569251779840e7',
        'week': 604800,
        'year': 31536000},

'unit_multipliers':
        {'bakers_dozen': 13,
        'dozen': 12,
        'gross': 144,
        'percent': '1/100'},

'velocity':
        {'knot': '463/900'},

'viscosity_absolute':
        {'poise': QQ((1, 10)),
        'reyn': '8896443230521/1290320000'},

'viscosity_kinematic':
        {'stokes': '1/10000'},

'viscosity_other':
        {'rhes': 10},

'volume':
        {'bag': '660732565629/6250000000000',
        'barrel': '9936705933/62500000000',
        'board_foot': '18435447/7812500000',
        'bucket': '473176473/31250000000',
        'bushel': '220244188543/6250000000000',
        'butt': '29810117799/62500000000',
        'cord': '884901456/244140625',
        'cubic_meter': 1,
        'cup': '473176473/2000000000000',
        'ephah': '1982197696887/50000000000000',
        'fifth': '473176473/625000000000',
        'firkin': '4091481/100000000',
        'fluid_dram': '473176473/128000000000000',
        'fluid_ounce': '473176473/16000000000000',
        'gallon': '473176473/125000000000',
        'gill': '473176473/4000000000000',
        'hogshead': '29810117799/125000000000',
        'imperial_gallon': '454609/100000000',
        'imperial_pint': '454609/800000000',
        'jeroboam': '473176473/156250000000',
        'jigger': '1419529419/32000000000000',
        'liter': '1/1000',
        'magnum': '473176473/250000000000',
        'minim': '157725491/2560000000000000',
        'noggin': '473176473/4000000000000',
        'omer': '1982197696887/500000000000000',
        'peck': '220244188543/25000000000000',
        'pint': '473176473/1000000000000',
        'pony': '1419529419/64000000000000',
        'puncheon': '9936705933/31250000000',
        'quart': '473176473/500000000000',
        'register_ton': '55306341/19531250',
        'seam': '220244188543/781250000000',
        'shot': '473176473/16000000000000',
        'stere': 1,
        'tablespoon': '473176473/32000000000000',
        'teaspoon': '157725491/32000000000000',
        'tun': '29810117799/31250000000',
        'uk_gallon': '454609/100000000',
        'uk_pint': '454609/800000000',
        'wine_bottle': '3/4000'}
}

unit_to_type = {}
value_to_unit = {}


def evalunitdict():
    """
    Replace all the string values of the unitdict variable by their
    evaluated forms, and builds some other tables for ease of use.
    This function is mainly used internally, for efficiency (and
    flexibility) purposes, making it easier to describe the units.

    EXAMPLES::

        sage: sage.symbolic.units.evalunitdict()
    """
    from sage.misc.sage_eval import sage_eval
    for key, value in unitdict.items():
        unitdict[key] = {a: (sage_eval(repr(b)) if isinstance(b, str) else b)
                         for a, b in value.items()}

    # FEATURE IDEA: create a function that would allow users to add
    # new entries to the table without having to know anything about
    # how the table is stored internally.

    #
    # Format the table for easier use.
    #
    for k, v in unitdict.items():
        for a in v:
            unit_to_type[a] = k

    for w in unitdict:
        for j in unitdict[w]:
            if isinstance(unitdict[w][j], tuple):
                unitdict[w][j] = unitdict[w][j][0]
        value_to_unit[w] = {b: a for a, b in unitdict[w].items()}


###############################################################################
#  Documentation for individual units.
#  Appears in unit's docstring.
###############################################################################

unit_docs = {
'acceleration_docs':
        {'gal': 'Abbreviation for galileo.\nDefined to be 1/100 meter/second^2.',
        'galileo': 'Defined to be 1/100 meter/second^2.',
        'gravity': 'Also called standard gravity.\nPhysical constant defined to be 9.80665 meter/second^2.'},

'amount_of_substance_docs':
        {'elementary_entity': 'Defined to be one elementary unit of choice, usually atoms or other elementary particles.\nApproximately equal to 1.6605e-24 moles.',
        'mole': 'SI base unit of quantity.\nDefined to be the amount of substance that has an equal number of elementary entities as there are atoms in 12 grams of carbon-12.\nEquivalent to Avogadros constant elementary entities or approximately equal to 6.022*10^23 elementary entities.'},

'angles_docs':
        {'arc_minute': 'Defined to be 1/60 of a degree or pi/10800 radians.',
        'arc_second': 'Defined to be 1/3600 of a degree or pi/648000 radians.',
        'degree': 'Defined to be pi/180 radians.',
        'grade': 'Defined to be pi/200 radians.',
        'quadrant': 'Equivalent to a right angle.\nDefined to be pi/2 radians.',
        'radian': 'SI derived unit of angle.\nDefined to be the angle subtended at the center of a circle by an arc that is equal in length to the radius of the circle.',
        'right_angle': 'Equivalent to a quadrant.\nDefined to be pi/2 radians.'},

'area_docs':
        {'acre': 'Defined to be 10 square chains or 4840 square yards.\nApproximately equal to 4046.856 square meters.',
        'are': 'Defined to be 100 square meters.',
        'barn': 'Defined to be 100 square femtometers or 10^-28 square meters.',
        'hectare': 'Defined to be 10000 square meters.',
        'rood': 'Defined to be 1/4 of an acre.\nApproximately equal to 1011.714 square meters.',
        'section': 'Equivalent to a square mile.\nApproximately equal to 2.59*10^6 square meters.',
        'square_chain': 'Defined to be 4356 square feet.\nApproximately equal to 404.9856 square meters.',
        'square_meter': 'SI derived unit of area.\nDefined to be meter^2.',
        'township': 'Defined to be 36 square miles.\nApproximately equal to 9.324*10^7 square meters.'},

'capacitance_docs':
        {'abfarad': 'Defined to be 10^9 farads.',
        'farad': 'SI derived unit of capacitance.\nDefined to be the charge in coulombs a capacitor will accept for the potential across it to change one volt.\nEquivalent to coulomb/volt.',
        'statfarad': 'CGS unit defined to be statcoulomb/statvolt.\nApproximately equal to 1.11265*10^-12 farads.'},

'charge_docs':
        {'abcoulomb': 'CGS unit defined to be 10 coulombs.',
        'coulomb': 'SI derived unit of charge.\nDefined to be the amount of electric charge transported by 1 ampere in 1 second.',
        'elementary_charge': 'Defined to be the amount of electric charge carried by a single proton or negative charge carried by a single electron.\nApproximately equal to 1.602176462*10^-19 coulombs.',
        'faraday': 'Defined to be the magnitude of electric charge in one mole of electrons.\nApproximately equal to 96485.3399 coulombs.',
        'franklin': 'CGS unit defined to be the amount of electric charge necessary such that if two stationary objects placed one centimeter apart had one franklin of charge each they would repel each other with a force of one dyne.\nApproximately equal to 3.3356*10^-10 coulombs.',
        'statcoulomb': 'Equivalent to franklin.\nApproximately equal to 3.3356*10^-10 coulombs.'},

'conductance_docs':
        {'abmho': 'Defined to be 10^9 siemens.',
        'mho': 'Equivalent to siemens.',
        'siemens': 'SI derived unit of conductance.\nDefined to be an ampere per volt or 1/ohm.'},

'current_docs':
        {'abampere': 'CGS unit defined to be 10 amperes.',
        'amp': 'Abbreviation for ampere.',
        'ampere': 'SI base unit of current.\nDefined to be the constant current which will produce an attractive force of 2*10^-7 newtons per meter between two straight, parallel conductors of infinite length and negligible circular cross section placed one meter apart in free space.',
        'biot': 'Equivalent to abampere.\nEqual to 10 amperes.',
        'statampere': 'CGS unit defined to be statcoulomb/second.\nApproximately equal to 3.335641*10^-10 amperes.'},

'electric_potential_docs':
        {'abvolt': 'Defined to be 10^-8 volts.',
        'statvolt': 'CGS unit defined to be the speed of light in a vacuum/10^6 volts or approximately 299.792 volts.',
        'volt': 'SI derived unit of electric potential.\nDefined to be the value of voltage across a conductor when a current of one ampere dissipates one watt of power.'},

'energy_docs':
        {'british_thermal_unit': 'Defined to be the amount of energy required to raise the temperature of one pound of liquid water from 60 degrees Fahrenheit to 61 degrees Fahrenheit at a constant pressure of one atmosphere.\nApproximately equal to 1055.05585 joules.',
        'btu': 'Abbreviation for British thermal unit.\nApproximately equal to 1055.05585 joules.',
        'calorie': 'Defined to be the amount of energy required to raise the temperature of one gram of liquid water one degree Celsius.\nEqual to 4.1868 joules.',
        'electron_volt': 'Defined to be the amount of kinetic energy gained by a single unbound electron when it accelerates through an electrostatic potential difference of 1 volt.\nApproximately equal to 1.602*10^-19 joules.',
        'erg': 'CGS unit for energy defined to be gram*centimeter^2/second^2.\nEqual to 10^-7 joules.',
        'ev': 'Abbreviation for electron volt.\nApproximately equal to 1.602*10^-19 joules.',
        'joule': 'SI derived unit of energy.\nDefined to be kilogram*meter^2/second^2.',
        'rydberg': 'Defined to be the absolute value of the binding energy of the electron in the ground state hydrogen atom.\nApproximately equal to 2.17987*10^-18 joules.',
        'therm': 'Defined to be 100,000 British thermal units.\nApproximately equal to 1.05505585*10^8 joules.'},

'fiber_linear_mass_density_docs':
        {'denier': 'Defined to be 1 gram per 9000 meters.\nEqual to 1/9000000 of a kilogram/meter.',
        'tex': 'Defined to be 1 gram per 1000 meters.\nEqual to 1/1000000 of a kilogram/meter.'},

'force_docs':
        {'dyne': 'CGS unit for force defined to be gram*centimeter/second^2.\nEqual to 10^-5 newtons.',
        'gram_weight': 'Defined to be the magnitude of the force exerted on one gram of mass by a 9.80665 meter/second^2 gravitational field.\nEqual to 1/1000 of a kilogram weight.\nEqual to 0.00980665 newtons.',
        'kilogram_force': 'Equivalent to a kilogram weight.\nEqual to 9.80665 newtons.',
        'kilogram_weight': 'Defined to be the magnitude of the force exerted on one kilogram of mass by a 9.80665 meter/second^2 gravitational field.\nEqual to 9.80665 newtons.',
        'newton': 'SI derived unit of force.\nDefined to be kilogram*meter/second^2.',
        'pound_force': 'Equivalent to a pound weight.\nApproximately equal to 4.44822 newtons.',
        'pound_weight': 'Defined to be the magnitude of the force exerted on one pound of mass by a 9.80665 meter/second^2 gravitational field.\nApproximately equal to 4.44822 newtons.',
        'poundal': 'Defined to be pound*foot/second^2.\nApproximately equal to 0.13825 newtons.',
        'ton_force': 'Defined to be 2000 pounds of force.\nApproximately equal to 8896.4432 newtons.'},

'frequency_docs':
        {'hertz': 'SI derived unit of frequency.\nDefined to be one complete cycle per second.'},

'illuminance_docs':
        {'foot_candle': 'Defined to be lumen/foot^2.\nApproximately equal to 10.764 lux.',
        'lux': 'SI derived unit of illuminance.\nDefined to be lumen/meter^2.',
        'phot': 'CGS unit defined to be 10000 lux.'},

'inductance_docs':
        {'abhenry': 'Defined to be 10^-9 henries.',
        'henry': 'SI derived unit of inductance./nDefined to be a volt per ampere per second.',
        'stathenry': 'CGS unit defined to be one statvolt*second/statampere.\nApproximately equal to 8.98758*10^11 henries.'},

'information_docs':
        {'bit': 'Base unit of information.\nDefined to be the maximum amount of information that can be stored by a device of other physical system that can normally exist in only two distinct states.',
        'byte': 'Defined to be 8 bits.',
        'nibble': 'Defined to be 4 bits.'},

'information_rate_docs':
        {'baud': 'Defined to be 1 bit/second.'},

'inverse_length_docs':
        {'diopter': 'Defined to be 1/meter.',
        'kayser': 'Defined to be 100/meter.'},

'length_docs':
        {'angstrom': 'Defined to be 10^-10 meters.',
        'astronomical_unit': 'Originally defined as the length of the semi-major axis of the elliptical orbit of the Earth around the Sun.\nRedefined for accuracy to be the radius of an unperturbed circular Newtonian orbit about the Sun of a particle having infinitesimal mass, moving with a mean motion of 0.01720209895 radians per day.\nApproximately equal to 1.496*10^11 meters.',
        'bolt': 'Defined to be 40 yards.\nEqual to 36.576 meters.',
        'cable_international': 'Nautical unit defined to be 1/10 of a nautical mile.\nEqual to 185.2 meters.',
        'cable_us': 'Nautical unit defined to be equal to 720 feet or 120 fathoms.\nEqual to 219.456 meters.',
        'caliber': 'Equal to 1/100 of an inch.\nEqual to 0.000254 meters.',
        'centimeter': 'Equal to 1/100 of a meter.',
        'chain': 'Surveying unit defined to be 66 feet.\nApproximately equal to 20.12 meters.',
        'cicero': 'Printing unit defined to be 12 didot points.\nApproximately equal to 0.004512 meters.',
        'cubit': 'Ancient unit of length defined to be 18 inches.\nEqual to 0.4572 meters.',
        'didot': 'Printing unit equal to 1/12 of a cicero.\nApproximately equal to 0.00037597 meters.',
        'dtp_point': 'The desktop publishing point is defined to be 1/72 of an inch.\nApproximately equal to 0.0003528 meters.',
        'ell': 'Ancient unit of length defined to be 45 inches.\nEqual to 1.143 meters.',
        'fathom': 'Nautical unit defined to be 6 feet.\nEqual to 1.8288 meters.',
        'feet': 'Equal to 12 inches.\nDefined to be 0.3048 meters.',
        'fermi': 'Equivalent to a femtometer.\nEqual to 10^-15 meters.',
        'foot': 'Equal to 12 inches.\nDefined to be 0.3048 meters.',
        'furlong': 'Defined to be 660 feet, or 1/8 of a mile.\nEqual to 201.168 meters.',
        'hand': 'Defined to be 4 inches.\nEqual to 0.1016 meters.',
        'inch': 'Equal to 1/12 of a foot.\nEqual to 0.0254 meters.',
        'kilometer': 'Equal to 1000 meters.\nEqual to 3280.8399 feet.',
        'league': 'Defined to be 3 miles.\nConventionally equal to the distance a person or horse can walk in one hour.\nEqual to 4828.032 meters.',
        'light_year': 'Defined to be the distance light travels in vacuum in 365.25 days.\nApproximately equal to 9.4607*10^15 meters.',
        'link': 'Surveying unit defined to be 1/100 of a chain.\nEqual to 0.201168 meters.',
        'meter': 'SI base unit of length.\nDefined to be the distance light travels in vacuum in 1/299792458 of a second.',
        'micron': 'Defined to be 10^-6 meters.',
        'mil': 'Defined to be 1/1000 of an inch.\nEqual to 0.0000254 meters.',
        'millimeter': 'Defined to be 1/1000 of a meter.\nEqual to 0.001 meters.',
        'mile': 'Defined to be 5280 feet.\nEqual to 1609.344 meters.',
        'nautical_mile': 'Nautical unit defined to be 1852 meters.',
        'parsec': 'Defined to be the length of the adjacent side of a right triangle whose angle is 1 arcsecond and opposite side equal to 1 astronomical unit, or 1 AU/arctan(1 arcsecond).\nApproximately equal to 30.857*10^15 meters.',
        'perch': 'Equivalent to rod.\nDefined to be 16.5 feet.\nEqual to 5.0292 meters.',
        'pica': 'Printing unit defined to be 12 dtp points.\nEqual to 1/72 of a foot.\nApproximately equal to 0.004233 meters.',
        'pole': 'Equivalent to rod.\nDefined to be 16.5 feet.\nEqual to 5.0292 meters.',
        'rod': 'Defined to be 16.5 feet.\nEqual to 5.0292 meters.',
        'rope': 'Defined to be 20 feet.\nEqual to 6.096 meters.',
        'skein': 'Defined to be 360 feet.\nEqual to 109.728 meters.',
        'stadion': 'Ancient unit of length defined to be 622 feet.\nEqual to 189.5856 meters.',
        'stadium': 'Defined to be 202 yards or 606 feet.\nEqual to 184.7088 meters.',
        'statute_mile': 'Equivalent to mile.\nDefined to be 5280 feet.\nEqual to 1609.344 meters.',
        'survey_foot': 'Defined to be 1200/3937 or approximately 0.3048006 meters.',
        'survey_mile': 'Defined to be 5280 survey feet.\nApproximately equal to 1609.347 meters.',
        'x_unit': 'Unit of length used to quote wavelengths of X-rays and gamma rays.\nApproximately equal to 1.0021*10^-13 meters.',
        'yard': 'Defined to be 3 feet.\nEqual to 0.9144 meters.'},

'luminance_docs':
        {'apostilb': 'Defined to be 10^-4 lamberts.\nEqual to 1/pi*candela/meter^2.',
        'lambert': 'Defined to be 10^4/pi candela/meter^2.',
        'nit': 'Equivalent to candela/meter^2.',
        'stilb': 'CGS unit equal to 10000 candela/meter^2.'},

'luminous_energy_docs':
        {'lumerg': 'Equivalent to lumen*second',
        'talbot': 'Equivalent to lumen*second.'},

'luminous_flux_docs':
        {'lumen': 'SI derived unit of luminous flux.\nDefined to be candela*steradian.'},

'luminous_intensity_docs':
        {'candela': 'SI base unit of luminous intensity.\nDefined to be the luminous intensity, in a given direction, of a source that emits monochromatic radiation of frequency 540*10^12 hertz and that has a radiant intensity in that direction of 1/683 watt per steradian.',
        'candle': 'Equivalent to candela.',
        'hefnerkerze': 'Old German unit defined to be a 8 millimeter wick burning amyl acetate with a flame height of 40 millimeters.\nApproximately equal to 0.9034 candelas.'},

'magnetic_field_docs':
        {'gauss': 'CGS unit defined to be a maxwell/centimeter^2.\nEqual to 1/10000 of a  tesla.',
        'tesla': 'SI derived unit of magnetic field.\nDefined to be the magnitude of a magnetic field such that a particle with a charge of 1 coulomb passing through that field at 1 meter/second will experience a force of 1 newton.'},

'magnetic_flux_docs':
        {'maxwell': 'CGS unit defined to be a gauss*centimeter^2 or 10^-8 webers.',
        'weber': 'SI derived unit of magnetic flux.\nDefined to be a change in magnetic flux of 1 weber per second will induce an electromotive force of 1 volt.'},

'magnetic_intensity_docs':
        {'oersted': 'CGS unit defined to be 1000/(4*pi) amperes per meter of flux path.'},

'magnetic_moment_docs':
        {'bohr_magneton': 'Physical constant defined to be the magnetic moment of an electron, or elementary_charge*h_bar/2*electron_rest_mass.\nApproximately equal to 9.274*10^-24 joules/tesla.',
        'nuclear_magneton': 'Physical constant defined to be the magnetic moment of a proton, or elementary_charge*h_bar/2*proton_rest_mass.\nApproximately equal to 5.05078324*10^-27 joules/tesla.'},

'magnetomotive_force_docs':
        {'ampere_turn': 'SI derived unit of magnetomotive force.\nDefined to be a direct current of 1 ampere flowing through a single turn loop in a vacuum.',
        'gilbert': 'CGS unit defined to be 10/(4*pi) ampere turns.'},

'mass_docs':
        {'amu': 'Abbreviation for atomic mass unit.\nApproximately equal to 1.660538782*10^-27 kilograms.',
        'assay_ton': 'Defined to be milligram*short_ton/ounce_troy.\nEqual to 7/240 of a kilogram.',
        'atomic_mass_unit': 'Defined to be one twelfth of the mass of an isolated atom of carbon-12 at rest and in its ground state.\nApproximately equal to 1.660538782*10^-27 kilograms.',
        'avoirdupois_ounce': 'Equivalent to ounce.\nEqual to 1/16 of an avoirdupois pound.\nApproximately equal to 0.02835 kilograms.',
        'avoirdupois_pound': 'Equivalent to pound.\nEqual to 16 avoirdupois ounces.\nApproximately equal to 0.45359 kilograms.',
        'bale': 'Equal to 500 pounds.\nApproximately equal to 226.796 kilograms.',
        'carat': 'Defined to be equal to 200 milligrams.\nCommonly denoted ct.',
        'cental': 'Equal to 100 pounds.\nApproximately equal to 45.36 kilograms.',
        'dalton': 'Equivalent to atomic_mass_unit.\nApproximately equal to 1.660538782*10^-27 kilograms.',
        'drachma': 'Ancient Greek unit of mass.\nEqual to 6 obols.\nApproximately equal to 0.00429234 kilograms.',
        'geepound': 'Equivalent to slug.\nApproximately equal to 14.5939 kilograms.',
        'grain': 'Historically based on the average mass of a single seed of a typical cereal.\nDefined in 1958 to be 64.79891 milligrams.',
        'gram': 'Equal to 0.0001 kilograms.',
        'gross_hundredweight': 'Equivalent to hundredweight.\nEqual to 112 pounds.\nApproximately equal to 50.802 kilograms.',
        'hundredweight': 'Defined to be 112 pounds.\nApproximately equal to 50.802 kilograms.',
        'kilogram': 'SI base unit of mass.\nDefined to be equal to the mass of the International Prototype Kilogram.\nAlmost exactly equal to the amount of mass in one liter of water.',
        'libra': 'Ancient Roman unit of mass.\nApproximately equal to 0.325971 kilogram.',
        'long_ton': 'Defined to be 2240 pounds.\nApproximately equal to 1016.05 kilograms.',
        'metric_ton': 'Defined to be 1000 kilograms.',
        'mina': 'Ancient Greek unit of mass.\nEqual to 100 drachma.\nApproximately equal to 0.429234 kilograms.',
        'net_hundredweight': 'Equivalent to cental.\nEqual to 100 pounds.\nApproximately equal to 45.36 kilograms.',
        'obol': 'Ancient Greek unit of mass.\nEqual to 1/6 of drachma.\nApproximately equal to 0.00071538 kilograms.',
        'ounce': 'Equal to 1/16 of pound.\nCommonly abbreviated oz.\nApproximately equal to 0.02835 kilograms.',
        'ounce_troy': 'Equal to 1/12 of pound_troy.\nApproximately equal to 0.031103 kilograms.',
        'pennyweight': 'Equal to 1/20 of ounce_troy.\nCommonly abbreviated dwt.\nApproximately equal to 0.001555 kilograms.',
        'pondus': 'Ancient Roman unit of mass.\nApproximately equal to 0.325969 kilograms.',
        'pound': 'Equal to 16 ounces.\nDefined to be exactly 0.45359237 kilograms.',
        'pound_troy': 'Equal to 12 ounce_troy.\nApproximately equal to 0.37324 kilograms.',
        'quintal': 'Equal to 100 kilograms.',
        'shekel': 'Ancient Hebrew unit of mass.\nApproximately equal to 0.0141 kilograms.',
        'short_hundredweight': 'Equivalent to cental.\nEqual to 100 pounds.\nApproximately equal to 45.36 kilograms.',
        'short_ton': 'Equivalent to ton.\nEqual to 2000 pounds.\nApproximately equal to 907.18 kilograms.',
        'slug': 'Defined to be a mass that is accelerated 1 ft/s^2 when 1 pound_force is exerted on it.\nApproximately equal to 14.5939 kilograms.',
        'solar_mass': 'Defined to be the mass of the Sun.\nAbout 332,950 times the mass of the Earth or 1,048 times the mass of Jupiter.\nApproximately equal to 1.98892*10^30 kilograms.',
        'stone': 'Defined to be 14 pounds.\nApproximately equal to 6.35 kilograms.',
        'talent': 'Ancient Greek unit of mass.\nEqual to 6000 drachmae.\nApproximately equal to 25.754 kilograms.',
        'ton': 'Equal to 2000 pounds.\nApproximately equal to 907.18 kilograms.',
        'tonne': 'Equivalent to metric_ton.\nDefined to be 1000 kilograms.',
        'wey': 'Defined to be 252 pounds.\nApproximately equal to 114.305 kilograms.'},

'power_docs':
        {'cheval_vapeur': 'Defined to be 75 kilogram force*meter/second.\nAlso known as metric horsepower.\nEqual to 735.49875 watts.',
        'horsepower': 'Defined to be 550 feet*pound force/second.\nApproximately equal to 745.7 watts.',
        'watt': 'SI derived unit of power.\nDefined to be joule/second or, in base units, kilogram*meter^2/second^3.'},

'pressure_docs':
        {'atmosphere': 'Defined to be 101325 pascals.',
        'bar': 'Defined to be 100000 pascals.',
        'barye': 'CGS unit defined to be dyne/centimeter^2.\nEqual to 1/10 of a pascal.',
        'inch_mercury': 'Defined to be 13595.1 kilogram/meter^3*inch*gravity.\nApproximately equal to 3386.389 pascals.',
        'millimeter_mercury': 'Defined to be 13595.1 kilogram/meter^3*millimeter*gravity.\nApproximately equal to 133.3224 pascals.',
        'mmhg': 'Abbreviation for millimeter mercury.\nApproximately equal to 133.3224 pascals.',
        'pa': 'Abbreviation for pascal.',
        'pascal': 'SI derived unit of pressure.\nDefined to be newton/meter^2 or, in base units, kilogram/(meter*second^2).',
        'pounds_per_square_inch': 'Defined to be pound force/inch^2.\nApproximately equal to 6894.76 pascals.',
        'psi': 'Abbreviation for pounds per square inch.\nApproximately equal to 6894.76 pascals.',
        'torr': 'Defined to be 1/760 of an atmosphere.\nApproximately equal to 133.322 pascals.'},

'radiation_absorbed_docs':
        {'gray': 'SI derived unit of absorbed radiation.\nDefined to be the absorption of one joule of ionizing radiation by one kilogram of matter.',
        'rad': 'Defined to be 1/100 of a gray.'},

'radiation_docs':
        {'becquerel': 'SI derived unit of radiation.\nDefined to be the activity of a quantity of radioactive material in which one nucleus decays per second.',
        'curie': 'Defined to be 37*10^9 becquerels.',
        'rutherford': 'Defined to be 10^6 becquerels.'},

'radiation_ionizing_docs':
        {'roentgen': 'Defined to be .000258 coulombs/kilogram.',
        'rontgen': 'Equivalent to roentgen.\nDefined to be .000258 coulombs/kilogram.'},

'resistance_docs':
        {'abohm': 'Defined to be 10^-9 ohms.',
        'ohm': 'SI derived unit of resistance.\nDefined to be a volt per ampere.',
        'statohm': 'CGS unit defined to be statvolt/statampere.\nApproximately equal to 8.98758*10^11 ohms.'},

'solid_angle_docs':
        {'steradian': 'SI derived unit of solid angle.\nDefined to be the solid angle subtended at the center of a sphere of radius r by a portion of the surface of the sphere having an area of r^2.'},

'temperature_docs':
        {'celsius': 'Defined to be -273.15 at absolute zero and 0.01 at the triple point of Vienna Standard Mean Ocean Water.\nCelsius is related to kelvin by the equation K = 273.15 + degrees Celsius.\nA change of 1 degree Celsius is equivalent to a change of 1 degree kelvin.',
        'centigrade': 'Equivalent to celsius.',
        'fahrenheit': 'Defined to be 32 degrees at the freezing point of water and 212 degrees at the boiling point of water, both at standard pressure (1 atmosphere).\nFahrenheit is related to kelvin by the equation K = 5/9*(degrees Fahrenheit + 459.67).\nA change of 1 degree fahrenheit is equal to a change of 5/9 kelvin.',
        'kelvin': 'SI base unit of temperature.\nDefined to be exactly 0 at absolute zero and 273.16 at the triple point of Vienna Standard Mean Ocean Water.',
        'rankine': 'Defined to be 0 at absolute zero and to have the same degree increment as Fahrenheit.\nRankine is related to kelvin by the equation K = 5/9*R.'},

'time_docs':
        {'century': 'Defined to be 100 years.\nEqual to 3153600000 seconds.',
        'day': 'Defined to be 24 hours.\nEqual to 86400 seconds.',
        'decade': 'Defined to be 10 years.\nEqual to 315360000 seconds.',
        'fortnight': 'Defined to be 2 weeks or 14 days.\nEqual to 1209600 seconds.',
        'hour': 'Defined to be 60 minutes.\nEqual to 3600 seconds.',
        'millenium': 'Defined to be 1000 years.\nEqual to 31536000000 seconds.',
        'minute': 'Defined to be 60 seconds.',
        'month': 'Defined to be 30 days.\nEqual to 2628000 seconds.',
        'second': 'SI base unit of time.\nDefined to be the duration of 9,192,631,770 periods of the radiation corresponding to the transition between the two hyperfine levels of the ground state of the caesium 133 atom.',
        'sidereal_day': 'Defined to be the time it takes for the Earth to make one complete rotation relative to the stars.\nApproximately equal to 86164.09 seconds.',
        'sidereal_second': 'Defined to be 1/86400 of a sidereal day.\nApproximately equal to 0.997269566329086 seconds.',
        'sidereal_year': 'Defined to be the time taken by the Earth to orbit the Sun once with respect to the fixed stars.\nApproximately equal to 31558149.7632 seconds.',
        'tropical_year': 'Defined to be the length of time that the Sun takes to return to the same position in the cycle of seasons, as seen from the Earth.\nApproximately equal to 31556925.1779840 seconds.',
        'week': 'Defined to be 7 days.\nEqual to 604800 seconds.',
        'year': 'Defined to be 365 days.\nEqual to 31536000 seconds.'},

'unit_multipliers_docs':
        {'bakers_dozen': 'Defined to be 13 items.',
        'dozen': 'Defined to be 12 items.',
        'gross': 'Defined to be 144 items.',
        'percent': 'Defined to be 1/100 of a quantity.'},

'velocity_docs':
        {'knot': 'Nautical unit of velocity defined to be a nautical mile per hour.\nApproximately equal to 0.5144 meter/second.'},

'viscosity_absolute_docs':
        {'poise': 'CGS unit defined to be 1/10 of pascal*second.',
        'reyn': 'Defined to be a pound_force*second/inch^2.\nApproximately equal to 6894.76 pascal*second.'},

'viscosity_kinematic_docs':
        {'stokes': 'CGS unit defined to be 1/10000 of meter^2/second.'},

'viscosity_other_docs':
        {'rhes': 'Defined to be 1/poise or 10/(pascal*second).'},

'volume_docs':
        {'bag': 'Defined to be 3 bushels.\nApproximately equal to 0.10572 cubic meters.',
        'barrel': 'Defined to be 42 gallons.\nApproximately equal to 0.15899 cubic meters.',
        'board_foot': 'Defined to be 144 cubic inches.\nApproximately equal to 0.0023597 cubic meters.',
        'bucket': 'Defined to be 4 gallons.\nApproximately equal to 0.0151416 cubic meters.',
        'bushel': 'Defined to be 2150.42 cubic inches.\nEquivalent to 4 pecks.\nApproximately equal to 0.035239 cubic meters.',
        'butt': 'Old English unit of wine casks defined to be 2 hogsheads or 126 gallons.\nApproximately equal to 0.476962 cubic meters.',
        'cord': 'Defined to be 8 feet x 8 feet x 4 feet.\nApproximately equal to 3.624556 cubic meters.',
        'cubic_meter': 'SI derived unit of volume.\nDefined to be meter^3.',
        'cup': 'Defined to be 8 fluid ounces.\nApproximately equal to 0.000236588 cubic meters.',
        'ephah': 'Ancient Hebrew unit of volume equal to 10 omers.\nApproximately equal to 0.03964 cubic meters.',
        'fifth': 'Defined to be 1/5 of a gallon.\nApproximately equal to 0.00075708 cubic meters.',
        'firkin': 'Defined to be 9 imperial gallons.\nApproximately equal to 0.04091 cubic meters.',
        'fluid_dram': 'Defined to be 1/8 of a fluid ounce.\nApproximately equal to 3.69669*10^-6 cubic meters.',
        'fluid_ounce': 'Defined to be 1/128 of a gallon.\nApproximately equal to 0.000029574 cubic meters.',
        'gallon': 'Defined to be 231 cubic inches.\nApproximately equal to 0.0037854 cubic meters.',
        'gill': 'Defined to be 4 fluid ounces.\nApproximately equal to 0.00011829 cubic meters.',
        'hogshead': 'Old English unit of wine casks defined to be 63 gallons.\nApproximately equal to 0.23848 cubic meters.',
        'imperial_gallon': 'Defined to be 4.54609 liters.\nEqual to 0.00454609 cubic meters.',
        'imperial_pint': 'Defined to be 1/8 of an imperial gallon.\nApproximately equal to 0.00056826 cubic meters.',
        'jeroboam': 'Defined to be 4/5 of a gallon.\nApproximately equal to 0.0030283 cubic meters.',
        'jigger': 'Defined to be 1 1/2 fluid ounces.\nApproximately equal to 0.00004436 cubic meters.',
        'liter': 'Defined to be 1 decimeter^3.\nEqual to 1/1000 of a cubic meter.',
        'magnum': 'Defined to be 1/2 a gallon.\nApproximately equal to 0.0018927 cubic meters.',
        'minim': 'Defined to be 1/480 of a fluid ounce.\nApproximately equal to 6.16115*10^-8 cubic meters.',
        'noggin': 'Equivalent to gill.\nDefined to be 4 fluid ounces.\nApproximately equal to 0.00011829 cubic meters.',
        'omer': 'Ancient Hebrew unit of volume equal to 9/20 of a peck.\nApproximately equal to 0.0039644 cubic meters.',
        'peck': 'Defined to be 1/4 of a bushel.\nApproximately equal to 0.0088098 cubic meters.',
        'pint': 'Defined to be 1/8 of a gallon.\nApproximately equal to 0.00047318 cubic meters.',
        'pony': 'Defined to be 3/4 of a fluid ounce.\nApproximately equal to 0.00002218 cubic meters.',
        'puncheon': 'Old English unit of wine casks defined to be 84 gallons.\nApproximately equal to 0.31797 cubic meters.',
        'quart': 'Defined to be 1/4 of a gallon.\nApproximately equal to 0.00094635 cubic meters.',
        'register_ton': 'Defined to be 100 cubic feet.\nApproximately equal to 2.83168 cubic meters.',
        'seam': 'Defined to be 8 bushels.\nApproximately equal to 0.281913 cubic meters.',
        'shot': 'Defined to be 1 fluid ounce.\nApproximately equal to 0.000029574 cubic meters.',
        'stere': 'Equivalent to cubic meter.',
        'tablespoon': 'Defined to be 1/2 of a fluid ounce.\nApproximately equal to 0.000014787 cubic meters.',
        'teaspoon': 'Defined to be 1/6 of a fluid ounce.\nEqual to 1/3 of a tablespoon.\nApproximately equal to 4.9289*10^-6 cubic meters.',
        'tun': 'Old English unit of wine casks defined to be 252 gallons.\nApproximately equal to 0.95392 cubic meters.',
        'uk_gallon': 'Equivalent to an imperial gallon.\nEqual to 0.00454609 cubic meters.',
        'uk_pint': 'Equivalent to and imperial pint.\nApproximately equal to 0.00056826 cubic meters.',
        'wine_bottle': 'Defined to be 750 milliliters.\nEqual to 0.00075 cubic meters.'}
}


###############################################################################
# Dictionary for converting from derived units to base SI units.
###############################################################################

unit_derivations = {'acceleration': 'length/time^2',
                    'area': 'length^2',
                    'capacitance': 'time^4*current^2/(length^2*mass)',
                    'charge': 'current*time',
                    'conductance': 'current^2*time^3/(mass*length^2)',
                    'electric_potential': 'mass*length^2/(current*time^3)',
                    'energy': 'mass*length^2/time^2',
                    'fiber_linear_mass_density': 'mass/length',
                    'force': 'mass*length/time^2',
                    'frequency': '1/time',
                    'illuminance': 'luminous_intensity*solid_angle/length^2',
                    'inductance': 'length^2*mass/(time^2*current^2)',
                    'information_rate': 'information/time',
                    'inverse_length': '1/length',
                    'luminance': 'luminous_intensity/length^2',
                    'luminous_energy': 'luminous_intensity*solid_angle*time',
                    'luminous_flux': 'luminous_intensity*solid_angle',
                    'magnetic_field': 'mass/(current*time^2)',
                    'magnetic_flux': 'mass*length^2/(current*time^2)',
                    'magnetic_intensity': 'current/length',
                    'magnetic_moment': 'current*length^2',
                    'power': 'mass*length^2/time^3',
                    'pressure': 'mass/(length*time^2)',
                    'radiation': '1/time',
                    'radiation_absorbed': 'length^2/time^2',
                    'radiation_ionizing': 'current*time/mass',
                    'resistance': 'mass*length^2/(current^2*time^3)',
                    'velocity': 'length/time',
                    'viscosity_absolute': 'mass/(length*time)',
                    'viscosity_kinematic': 'length^2/time',
                    'viscosity_other': 'length*time/mass',
                    'volume': 'length^3'
                    }


def vars_in_str(s):
    """
    Given a string like ``'mass/(length*time)'``, return the list
    ``['mass', 'length', 'time']``.

    INPUT:

    - ``s`` -- string

    OUTPUT: list of strings (unit names)

    EXAMPLES::

        sage: sage.symbolic.units.vars_in_str('mass/(length*time)')
        ['mass', 'length', 'time']
    """
    return re.findall('[a-z|_]+', s)


def unit_derivations_expr(v):
    """
    Given derived units name, returns the corresponding units expression.

    For example, given ``'acceleration'`` output the symbolic expression
    ``length/time^2``.

    INPUT:

    - ``v`` -- string; name of a unit type such as ``'area'``, ``'volume'``, etc.

    OUTPUT: a symbolic expression

    EXAMPLES::

        sage: sage.symbolic.units.unit_derivations_expr('volume')
        length^3
        sage: sage.symbolic.units.unit_derivations_expr('electric_potential')
        length^2*mass/(current*time^3)

    If the unit name is unknown, a :exc:`KeyError` is raised::

        sage: sage.symbolic.units.unit_derivations_expr('invalid')
        Traceback (most recent call last):
        ...
        KeyError: 'invalid'
    """
    v = str(v)
    Z = unit_derivations[v]
    if isinstance(Z, str):
        d = {x: str_to_unit(x) for x in vars_in_str(Z)}
        from sage.misc.sage_eval import sage_eval
        Z = sage_eval(Z, d)
        unit_derivations[v] = Z
    return Z


@instancedoc
class UnitExpression(Expression):
    """
    A symbolic unit.

    EXAMPLES::

        sage: acre = units.area.acre
        sage: type(acre)
        <class 'sage.symbolic.units.UnitExpression'>

    TESTS::

        sage: bool(loads(dumps(acre)) == acre)
        True
        sage: type(loads(dumps(acre)))
        <class 'sage.symbolic.units.UnitExpression'>
    """
    def _instancedoc_(self):
        """
        Return docstring for this unit.

        EXAMPLES::

            sage: print(units.area.acre.__doc__)
            Defined to be 10 square chains or 4840 square yards.
            Approximately equal to 4046.856 square meters.
        """
        return unitdocs(self)


def str_to_unit(name):
    """
    Create the symbolic unit with given name.  A symbolic unit is a
    class that derives from symbolic expression, and has a specialized
    docstring.

    INPUT:

    - ``name`` -- string

    OUTPUT: a :class:`UnitExpression`

    EXAMPLES::

        sage: sage.symbolic.units.str_to_unit('acre')
        acre
        sage: type(sage.symbolic.units.str_to_unit('acre'))
        <class 'sage.symbolic.units.UnitExpression'>
    """
    return UnitExpression(SR, SR.var(name))


class Units(ExtraTabCompletion):
    """
    A collection of units of some type.

    EXAMPLES::

        sage: units.power
        Collection of units of power: cheval_vapeur horsepower watt
    """
    def __init__(self, data, name=''):
        """
        EXAMPLES::

            sage: sage.symbolic.units.Units(sage.symbolic.units.unitdict, 'all units')
            Collection of units of all units: acceleration ... volume
        """
        self.__name = name
        self.__data = data
        self.__units = {}

    def __getstate__(self):
        """
        Used for pickling.   We throw away all cached information.

        EXAMPLES::

            sage: type(units.__getstate__()[0])
            <... 'str'>
            sage: type(units.__getstate__()[1])
            <... 'dict'>
            sage: loads(dumps(units)) == units
            True
            sage: loads(dumps(units.area)) == units.area
            True
            sage: bool(loads(dumps(units.area.acre)) == units.area.acre)
            True
        """
        return (self.__name, self.__data)

    def __setstate__(self, state):
        """
        Used for unpickling.  See __getstate__.

        EXAMPLES::

            sage: state = units.__getstate__()
            sage: units.__setstate__(state)
        """
        self.__name = state[0]
        self.__data = state[1]
        self.__units = {}

    def __eq__(self, other):
        """
        Compare two collections of units, or a collection of units
        with some other object.

        EXAMPLES::

            sage: units.length == 10
            False
            sage: units.length == units.length
            True
            sage: units.length == units.mass
            False
        """
        if not isinstance(other, Units):
            return False
        return (self.__name, self.__data) == (other.__name, other.__data)

    def __ne__(self, other):
        """
        Test for unequality.

        EXAMPLES::

            sage: units.length != 5
            True
            sage: units.length != units.length
            False
            sage: units.length != units.mass
            True
        """
        return not (self == other)

    def _tab_completion(self):
        """
        Return tab completions.

        This complements the usual content of :func:`dir`, with the
        list of the names of the unit collections (resp. units) for
        :obj:`units` (resp. its subcollections), in particular for tab
        completion purposes.

        .. SEEALSO:: :class:`ExtraTabCompletion`

        EXAMPLES::

            sage: units.area._tab_completion()
            ['acre', 'are', 'barn', 'hectare', 'rood', 'section', 'square_chain', 'square_meter', 'township']
            sage: units._tab_completion()
            ['acceleration', ..., 'volume']
            sage: units.force._tab_completion()
            ['dyne', ..., 'ton_force']

            sage: dir(units)
            ['_Units__data', ..., 'acceleration', ..., 'volume']
            sage: dir(units.force)
            ['_Units__data', ..., 'dyne', ..., 'ton_force']
        """
        return sorted([x for x in self.__data if '/' not in x])

    def __getattr__(self, name):
        """
        Return the unit with the given name.

        EXAMPLES::

            sage: units.area
            Collection of units of area: acre are barn hectare rood section square_chain square_meter township
            sage: units.area.barn
            barn

        Units are cached::

            sage: units.area.acre is units.area.acre
            True
        """
        if name in self.__units:
            return self.__units[name]
        if len(unit_to_type) == 0:
            evalunitdict()
        try:
            v = self.__data[name]
        except KeyError:
            raise AttributeError
        if isinstance(v, dict):
            U = Units(self.__data[name], name)
        else:
            U = str_to_unit(name)
        self.__units[name] = U
        return U

    def __repr__(self):
        """
        Return string representation of this collection of units.

        EXAMPLES::

            sage: units.__repr__()
            'Collection of units: acceleration ... volume'
            sage: units.area.__repr__()
            'Collection of units of area: acre are barn hectare rood section square_chain square_meter township'
        """
        name = ' of ' + self.__name if self.__name else ''
        return "Collection of units{0}: {1}".format(name, ' '.join(sorted([str(x) for x in self.__data])))


units = Units(unitdict, '')


def unitdocs(unit):
    r"""
    Return docstring for the given unit.

    INPUT:

    - ``unit`` -- a unit

    OUTPUT: string

    EXAMPLES::

        sage: sage.symbolic.units.unitdocs('meter')
        'SI base unit of length.\nDefined to be the distance light travels in vacuum in 1/299792458 of a second.'
        sage: sage.symbolic.units.unitdocs('amu')
        'Abbreviation for atomic mass unit.\nApproximately equal to 1.660538782*10^-27 kilograms.'

    Units not in the list unit_docs will raise a :exc:`ValueError`::

        sage: sage.symbolic.units.unitdocs('earth')
        Traceback (most recent call last):
        ...
        ValueError: no documentation exists for the unit earth
    """
    if is_unit(unit):
        return unit_docs[unit_to_type[str(unit)] + "_docs"][str(unit)]
    raise ValueError("no documentation exists for the unit %s" % unit)


def is_unit(s) -> bool:
    """
    Return a boolean when asked whether the input is in the list of units.

    INPUT:

    - ``s`` -- an object

    OUTPUT: boolean

    EXAMPLES::

        sage: sage.symbolic.units.is_unit(1)
        False
        sage: sage.symbolic.units.is_unit(units.length.meter)
        True

    The square of a unit is not a unit::

        sage: sage.symbolic.units.is_unit(units.length.meter^2)
        False

    You can also directly create units using var, though they won't have
    a nice docstring describing the unit::

        sage: sage.symbolic.units.is_unit(var('meter'))
        True
    """
    return str(s) in unit_to_type


def convert(expr, target):
    """
    Convert units between ``expr`` and ``target``. If ``target`` is ``None``
    then converts to SI base units.

    INPUT:

    - ``expr`` -- the symbolic expression converting from

    - ``target`` -- (default: ``None``) the symbolic expression converting to

    OUTPUT: a symbolic expression

    EXAMPLES::

        sage: sage.symbolic.units.convert(units.length.foot, None)
        381/1250*meter
        sage: sage.symbolic.units.convert(units.mass.kilogram, units.mass.pound)
        100000000/45359237*pound

    This raises :exc:`ValueError` if expr and target are not convertible::

        sage: sage.symbolic.units.convert(units.mass.kilogram, units.length.foot)
        Traceback (most recent call last):
        ...
        ValueError: Incompatible units
        sage: sage.symbolic.units.convert(units.length.meter^2, units.length.foot)
        Traceback (most recent call last):
        ...
        ValueError: Incompatible units

    Recognizes derived unit relationships to base units and other derived units::

        sage: sage.symbolic.units.convert(units.length.foot/units.time.second^2, units.acceleration.galileo)
        762/25*galileo
        sage: sage.symbolic.units.convert(units.mass.kilogram*units.length.meter/units.time.second^2, units.force.newton)
        newton
        sage: sage.symbolic.units.convert(units.length.foot^3, units.area.acre*units.length.inch)
        1/3630*(acre*inch)
        sage: sage.symbolic.units.convert(units.charge.coulomb, units.current.ampere*units.time.second)
        (ampere*second)
        sage: sage.symbolic.units.convert(units.pressure.pascal*units.si_prefixes.kilo, units.pressure.pounds_per_square_inch)
        1290320000000/8896443230521*pounds_per_square_inch

    For decimal answers multiply 1.0::

        sage: sage.symbolic.units.convert(units.pressure.pascal*units.si_prefixes.kilo, units.pressure.pounds_per_square_inch)*1.0
        0.145037737730209*pounds_per_square_inch

    You can also convert quantities of units::

        sage: sage.symbolic.units.convert(cos(50) * units.angles.radian, units.angles.degree)
        degree*(180*cos(50)/pi)
        sage: sage.symbolic.units.convert(cos(30) * units.angles.radian, units.angles.degree).polynomial(RR)
        8.83795706233228*degree
        sage: sage.symbolic.units.convert(50 * units.length.light_year / units.time.year, units.length.foot / units.time.second)
        6249954068750/127*(foot/second)

    Quantities may contain variables (not for temperature conversion, though)::

        sage: sage.symbolic.units.convert(50 * x * units.area.square_meter, units.area.acre)
        acre*(1953125/158080329*x)
    """
    base_target = target
    z = {}
    tz = {}

    for x in expr.variables():
        if is_unit(x):
            if unit_to_type[str(x)] == 'temperature':
                return convert_temperature(expr, target)
            else:
                z[x] = base_units(x)

    expr = expr.subs(z)

    if target is None:
        return expr
    else:
        for y in base_target.variables():
            if is_unit(y):
                tz[y] = base_units(y)
        base_target = base_target.subs(tz)
        coeff = (expr / base_target).expand()

        for variable in coeff.variables():
            if is_unit(str(variable)):
                raise ValueError("Incompatible units")

        return coeff.mul(target, hold=True)


def base_units(unit):
    """
    Convert unit to base SI units.

    INPUT:

    - ``unit`` -- a unit

    OUTPUT: a symbolic expression

    EXAMPLES::

        sage: sage.symbolic.units.base_units(units.length.foot)
        381/1250*meter

    If unit is already a base unit, it just returns that unit::

        sage: sage.symbolic.units.base_units(units.length.meter)
        meter

    Derived units get broken down into their base parts::

        sage: sage.symbolic.units.base_units(units.force.newton)
        kilogram*meter/second^2
        sage: sage.symbolic.units.base_units(units.volume.liter)
        1/1000*meter^3

    Returns variable if 'unit' is not a unit::

        sage: sage.symbolic.units.base_units(var('x'))
        x
    """
    from sage.misc.sage_eval import sage_eval
    if str(unit) not in unit_to_type:
        return unit
    elif unit_to_type[str(unit)] == 'si_prefixes' or unit_to_type[str(unit)] == 'unit_multipliers':
        number = unitdict[unit_to_type[str(unit)]][str(unit)]
        return (number if number in QQ else sage_eval(number))
    else:
        v = SR.var(unit_to_type[str(unit)])
        if str(v) in unit_derivations:
            base = unit_derivations_expr(v)
            for i in base.variables():
                base = base.subs({i: SR.var(value_to_unit[str(i)][1])})
            number = unitdict[str(v)][str(unit)]
            return base * (number if number in QQ else sage_eval(number))
        else:
            base = SR.var(value_to_unit[str(v)][1])
            number = unitdict[str(v)][str(unit)]
            return base * (number if number in QQ else sage_eval(number))


def convert_temperature(expr, target):
    """
    Function for converting between temperatures.

    INPUT:

    - ``expr`` -- a unit of temperature
    - ``target`` -- a units of temperature

    OUTPUT: a symbolic expression

    EXAMPLES::

        sage: t = 32*units.temperature.fahrenheit
        sage: t.convert(units.temperature.celsius)
        0
        sage: t.convert(units.temperature.kelvin)
        273.150000000000*kelvin

    If target is None then it defaults to kelvin::

        sage: t.convert()
        273.150000000000*kelvin

    This raises :exc:`ValueError` when either input is not a unit of temperature::

        sage: t.convert(units.length.foot)
        Traceback (most recent call last):
        ...
        ValueError: cannot convert
        sage: wrong = units.length.meter*units.temperature.fahrenheit
        sage: wrong.convert()
        Traceback (most recent call last):
        ...
        ValueError: cannot convert

    We directly call the convert_temperature function::

        sage: sage.symbolic.units.convert_temperature(37*units.temperature.celsius, units.temperature.fahrenheit)
        493/5*fahrenheit
        sage: 493/5.0
        98.6000000000000
    """
    if len(expr.variables()) != 1:
        raise ValueError("cannot convert")
    elif target is None or unit_to_type[str(target)] == 'temperature':
        from sage.misc.sage_eval import sage_eval
        expr_temp = expr.variables()[0]
        coeff = expr / expr_temp
        if target is not None:
            target_temp = target.variables()[0]
        a = sage_eval(unitdict['temperature'][str(expr_temp)],
                      locals={'x': coeff})
        if target is None or target_temp == units.temperature.kelvin:
            return a[0] * units.temperature.kelvin
        elif target_temp == units.temperature.celsius or target_temp == units.temperature.centigrade:
            return a[1] * target_temp
        elif target_temp == units.temperature.fahrenheit:
            return a[2] * units.temperature.fahrenheit
        elif target_temp == units.temperature.rankine:
            return a[3] * target_temp
    else:
        raise ValueError("cannot convert")
