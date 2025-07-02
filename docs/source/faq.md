General
=======

Q: How can I find all the classes available in MOOSE?
-----------------------------------------------------

### A: In Python, after importing moose, run the command `moose.le('/classes')`

Q: How can I see the documentation on a class?
----------------------------------------------

### A: You can run `help(moose.{classname}) ~moose.doc('{classname}')` in Python.

Here `{classname}` should be replaced by the name of the class you want
to know about.

For example, to learn more about the `Compartment` class you can run
either of the following code snippets:

``` {.python}
import moose
moose.doc('Compartment')
```

or

``` {.python}
import moose
help(moose.Compartment)
```

Q: How can I see the documentation of a specific field of a class?
------------------------------------------------------------------

### A: Run `moose.doc('{classname}.{fieldname}')`

Q: What unit system does MOOSE use?
-----------------------------------

### A: MOOSE is unit agnostic. If all the parameters in your model are in a consistent unit system, results will be in the same system. We recommend that you use SI units throughout.

Experimentalists use different units for different quantities for
convenience. These often get mixed up in computational models and can
produce unexpected results due to inconsistent units. The physiological
unit system uses $mV$ for membrane potential, $ms$ for time. Now, if you
express membrane capacitance in $pF$, your membrane resistance must be
in in $G\Omega$ to get consistent time course in your simulation. This
is because only then your membrane time constant
($\tau = R (G\Omega) \times C (pF)$) will be in $ms$. But then your
currents would be in $pA$ ($I = \frac{V (mV) }{ R (G\Omega)}$).

Another source of frustration is that dimensions of neurons are usually
expressed in $\mu m$. But conveniently, the specific capacitance of cell
membrane is close to $1 uF/cm^{2}$. So if you are trying to set the
absolute capacitance of a compartment based on this, you must convert
its length and diameter to $cm$ first and compute the surface area in
$cm^{2}$.

Q: How can I find all elements of a certain type in a model directly under a given element?
-------------------------------------------------------------------------------------------

### A: Use `wildcardFind` in the form `moose.wildcardFind('{root}/#[TYPE={typename}]')` where `root` is the element under which to search.

Example: You loaded the [Kholodenko 2000
model](https://moose.ncbs.res.in/readthedocs/user/py/tutorials/ChemicalOscillators.html#slow-feedback-oscillator)
as `'/Kholodenko'` , then the container for the chemicals and reactions
is `'/Kholodenko/kinetics/MAPK'` . You want to search for elements of
type `Pool` (which represents a pool of molecules or ions) directly
under `MAPK`.

``` {.python wrap="yes"}
In [19]: moose.wildcardFind('/Kholodenko/kinetics[0]/MAPK[0]/#[TYPE==Pool]')
Out[19]: [<moose.Reac id=521 dataIndex=0 path=/Kholodenko[0]/kinetics[0]/MAPK[0]/Neg_feedback[0]>]
```

This returns a list of elements of type `Pool` that are children of
`/Kholodenko/kinetics/MAPK` . But not its grand children or beyond.

Q: How can I find all elements of a certain type at all depths in a model?
--------------------------------------------------------------------------

### A: Use `wildcardFind` in the form `moose.wildcardFind('{root}/##[TYPE={typename}]')`

Example: In the [Kholodenko 2000
model](https://moose.ncbs.res.in/readthedocs/user/py/tutorials/ChemicalOscillators.html#slow-feedback-oscillator)
you want to search for elements of type `Pool` (which represents a pool
of molecules or ions) you loaded the model as `"/Kholodenko\"`, then

``` {.python}
In [19]: moose.wildcardFind('/Kholodenko/##[TYPE==Pool]')
Out[19]: [<moose.Reac id=521 dataIndex=0 path=/Kholodenko[0]/kinetics[0]/MAPK[0]/Neg_feedback[0]>]
```

The `##` here tells `wildcardFind` to do a recursive seach, i.e.,
starting at `{root}` , look through its children, grand children, all
the way down to the leaves of the element tree.

Q: How can I find an element with a particular name in a large model?
---------------------------------------------------------------------

### A: Use `wildcardFind` in the form `moose.wildcardFind('{root}/##[FIELD(name)={name}]')`

This always returns a list of elements. If the name is unique under
root, it will be a list with a single element. If there are multiple
elements with the same name, this will return a list containing all
those elements. Example: In the [Kholodenko 2000
model](https://moose.ncbs.res.in/readthedocs/user/py/tutorials/ChemicalOscillators.html#slow-feedback-oscillator)
you want to search for the element `"Neg_feedback"`, and you loaded the
model as `"/Kholodenko\"`, then

``` {.python}
In [19]: moose.wildcardFind('/Kholodenko/##[FIELD(name)==Neg_feedback]')
Out[19]: [<moose.Reac id=521 dataIndex=0 path=/Kholodenko[0]/kinetics[0]/MAPK[0]/Neg_feedback[0]>]
```

If you search for `"MAPK"` :

``` {.python}
In [21]: moose.wildcardFind('/Kholodenko/##/MAPK')
Out[21]:
[<moose.Neutral id=489 dataIndex=0 path=/Kholodenko[0]/kinetics[0]/MAPK[0]>,
<moose.Pool id=491 dataIndex=0 path=/Kholodenko[0]/kinetics[0]/MAPK[0]/MAPK[0]>,
<moose.Table2 id=552 dataIndex=0 path=/Kholodenko[0]/data[0]/MAPK[0]>]
```

Note that the recursive search wildcard (`##`) must be separated by `/`
(slash). `moose.wildcardFind('/Kholodenko/##MAPK')` returns an empty
list.

Q: When trying to plot data from moose `Table` objects, I see this error \"`ValueError: x and y must have same first dimension, but have shapes ...`\". What is the problem?
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

``` {.example}
  File "C:\moose-core\tests\core\debug_hhchanf2d.py", line 86, in test_vclamp
    axes[1].plot(t, gktab.vector, label=f'({vstep * 1e3} mV)')
  File "C:\miniforge3\envs\track\Lib\site-packages\matplotlib\axes\_axes.py", line 1721, in plot
    lines = [*self._get_lines(self, *args, data=data, **kwargs)]
            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\miniforge3\envs\track\Lib\site-packages\matplotlib\axes\_base.py", line 303, in __call__
    yield from self._plot_args(
            ^^^^^^^^^^^^^^^^
  File "C:\miniforge3\envs\track\Lib\site-packages\matplotlib\axes\_base.py", line 499, in _plot_args
    raise ValueError(f"x and y must have same first dimension, but "
ValueError: x and y must have same first dimension, but have shapes (1050001,) and (2100002,)
```

### A: This may happen due to multiple connections.

In this instance, you are creating the time vector based on the length
of one `Table`, and plotting another. One of the tables (here the second
one) was doubly connected, i.e., its `requestOut` message was passed in
two connect calls, like this:

``` {.python}
moose.connect(gktab, 'requestOut', channel1, 'Gk')
moose.connect(gktab, 'requestOut', channel2, 'Gk')
```

You should check the length of the `vector` attribute of a table against
the number of time points expected from total simulation time (simtime
in `moose.start(simtime)` ) and the `dt` of the table. If the length of
`vector` is a multiple of the expected number, then you have multiple
connection issue with that table.

Electrophysiology
=================

Q: Does setting `length` and `diameter` of a `Compartment` update its membrane resistance `Rm` and capacitance `Cm` fields?
---------------------------------------------------------------------------------------------------------------------------

### A: No.

`length` and `diameter` are utility fields that can be used for
visualization and format conversion. They are not used by MOOSE for
computing `Rm` and `Cm`. In some other tools the user specifies the
specific resistance and specific capacitance of the membrane along with
`length` and `diameter` of a compartment. And the software computes the
total mebrane resistance and capacitance based on the surface area.
However, in MOOSE `Rm` and `Cm` represent total membrane resistance and
capacitance respectively, and are set directly. If you have a
cylindrical compartment with given length and diameter, and the specific
resistance and capacitance are `RM` and `CM` respectively, then you
should set:

``` {.python}
surface_area = pi * compartment.length * compartment.diameter
compartment.Rm = RM / surface_area
compartment.Cm = CM * surface_area
```

where `pi` is imported from any of the modules like `math` or `numpy`,
or set directly to 3.14159... You must ensure unit consistency between
length and specific resistance/capacitance here.

Q: What is the use of the field `Ra` in a `Compartment`?
--------------------------------------------------------

### A: `Ra` represents total axial resistance (also called cytoplasmic resistance) of the compartment.

It is relevant only for multicompartmental models, and represents the
resistance faced by current flowing from one compartment to the next.

Q: What does `HH` in `HHChannel`, `HHGate`, etc. stand for?
-----------------------------------------------------------

### A: `HH` stands for Hodgkin and Huxley.

They figured out the dynamics of ion channels that generate action
potentials in nerve fibres, and these models use their formulation to
implement channel dynamics. See
[here](https://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model),
[here](https://physoc.onlinelibrary.wiley.com/pb-assets/assets/14697793/Companion_Guide_V2-1659940264.pdf)
and [here](https://doi.org/10.1113/jphysiol.2012.230458) to learn more
about this.

Q: How do I model Hodgkin-Huxley type ion channels?
---------------------------------------------------

### A: Create instances of the `HHChannel` class for fast computations, or `HHChannelF` class for more flexible computations.

The `HHChannel` class is designed for faster simulations, and avoids
explicitly evaluating Hodgkin-Huxley-type expressions at each simulation
step. Instead, it finds the `alpha(V)` and `beta(V)` for membrane
voltage `V` using lookup/interpolation tables.

The `HHChannelF` class is for convenience and accuracy, and allows you
to explicitly set the gate expressions as strings. However, it evaluates
the formula at each time step, which is slower than table lookup in
`HHChannel`.

Q: How do I set up the channel gating equations in `HHChannel`?
---------------------------------------------------------------

### A: The gating parameters in `HHChannel` are computed via `HHGate` objects. They use tables mapping voltage values to gating parameter values to find the gating parameters for a given voltage.

1.  Set `HHChannel.Xpower`, `HHChannel.Ypower` corresponding to the
    number of gating particles in the channel model.
    -   For example, for the Na+ channel conductance in the squid giant
        axon is calculated at voltage $V$ as
        $G_{Na}(V) = \bar{G}_{Na} \times m(V)^{3} \times h(V)$. So set
        `channel.Xpower= 3 and ~channel.Ypower=1`.
    -   When these powers are set to positive values, corresponding
        `HHGate` objects are created.
2.  Next setup the lookup tables for gating parameters based on their
    equations using one of the following methods:
    1.  Use `HHGate.setupAlpha()` (for alpha-beta form of the HH
        equations) or `HHGate.setupTau()` (for the tau-inf form of the
        HH equations) function.
    2.  (Introduced in 2025-03 in development branch) Assign
        mathematical expression strings to the fields `HHGate.alphaExpr`
        and `HHGate.betaExpr` for alpha-beta form, or `HHGate.tauExpr`
        and `HHGate.infExpr` for tau-inf form
    3.  Explicitly set the `HHGate.tableA` and `HHGate.tableB` fields
        with arrays computed in Python.

Q: How do I use `HHGate.setupAlpha()` or `HHGate.setupTau()` for setting up the lookup tables of an `HHGate` object?
--------------------------------------------------------------------------------------------------------------------

### A: Use `HHGate.setupAlpha(...)` for alpha-beta form of the HH equations or `HHGate.setupTau(...)` for the tau-inf form of the HH equations.

These functions compute the HH-gate expressions in the generic form:

$y(x) = \frac{A + Bx}{C + exp((x+D)/F)}$

Thus for each of `alpha` and `beta` (or `tau` and `inf` depeding on the
formulation) there are 5 constant coefficients: A...F, some of which may
be `0`. The `setupAlpha()` function (or `setupTau()`) takes a list of 13
parameters: 5 each for the constants in the expressions for `alpha` and
`beta` (or `tau` and `inf`), number of divisions in the interpolation
table (`divs`), the lower bound of the interpolation table (`min`), and
the upper bound of the interpolation table (`max`).

Below is an example of setting up Hodgkin-Huxley\'s Na channel using
`setupAlpha`. In their formulation

$\alpha_{m} = \frac{0.1 (25 - V)}{exp(\frac{25 - V}{10}) - 1}$

Thus $A = 0.1 * 25$ , $B = -0.1$ , $C = -1$ , $D = -25$ and $F = -10$ .
Similarly,

$\beta_{m} = 4 exp(\frac{-V}{18})$

so, $A = 4$, $B = 0$, $C = 0$, $D = 0$ and $F = 18$ .

We expect the membrane voltage to stay within -110 mV to 50 mV under
physiological conditions, and divide this range into in 3000 points for
lookup. We also want the lookup to use linear interpolation for voltage
values falling between two entries in the table.

``` {.python}
import moose

chan = moose.HHChannel('channel')
chan.Xpower = 3  # this will also initialize the HHGate element channel/gateX
chan.Ypower = 1  # this will also initialize the HHGate element channel/gateY
m_gate = moose.element(f'{chan.path}/gateX')
h_gate = moose.element(f'{chan.path}/gateY')
vmin = -110
vmax = 50
vdivs = 3000
m_gate.setupAlpha([
    0.1 * 25.0,                # A_A
    -0.1,                      # A_B
    -1.0,                      # A_C
    -25.0,                     # A_D
    -10.0,                     # A_F
    4.0,                       # B_A
    0.0,                       # B_B
    0.0,                       # B_C
    0.0,                       # B_D
    18.0,                      # B_F
    vdivs,
    vmin,
    vmax])
m_gate.useInterpolation = True   # use linear interpolation instead of direct lookup

# Similarly for h_gate ...
```

Q: How do I use explicit mathematical expressions for setting up the gating dynamics of `HHGate`?
-------------------------------------------------------------------------------------------------

### (Introduced in 2025-03 in development branch) Assign mathematical expression strings to the fields `HHGate.alphaExpr` and `HHGate.betaExpr` for alpha-beta form, or `HHGate.tauExpr` and `HHGate.infExpr` for tau-inf form.

In this case you must explicitly set the `min`, `max`, and `divs` fields
of each gate. The expressions should use
[exprtk](https://github.com/ArashPartow/exprtk) syntax. Do not forget to
convert the expressions to reflect the unit system of your entire model.
We recommend adhering to SI units, but the code sample below shows
original Hodgkin-Huxley formulation in physiological units. Note that
the `min` and `max` voltages are in `mV`.

``` {.python}
import moose

chan = moose.HHChannel('channel')
chan.Xpower = 3  # this will also initialize the HHGate element channel/gateX
chan.Ypower = 1  # this will also initialize the HHGate element channel/gateY
m_gate = moose.element(f'{chan.path}/gateX')
h_gate = moose.element(f'{chan.path}/gateY')
m_gate.alphaExpr = '0.1 * (25 - v)/(exp((25 - v)/10) - 1)'
m_gate.betaExpr = '4 * exp(-v/18)'
h_gate.alphaExpr = '0.07 * exp(-v/20)'
h_gate.betaExpr = '1/(exp((30-v)/10) + 1)'
for gate in (m_gate, h_gate):
   gate.useInterpolation = True
   gate.divs = 1000
   gate.min = -30.0
   gate.max = 120.0
```

Q: How do I explicitly set the lookup tables for `HHGates` when using an `HHChannel`?
-------------------------------------------------------------------------------------

### A: Compute the `alpha(V)` and `beta(V)` values for the desired range of voltages.

Then for the target gate, assign `alpha(V)` series to the `tableA`
field, and `alpha(V) + beta(V)` to the `tableB` field. If you have few
voltage values, or if you want higher accuracy, set `useInterpolation`
to `True`. You must explicitly set the `min`, `max`, and `divs` fields
of each gate. Note that we are using SI units here and the HH equations
have been modified accordingly.

``` {.python}
import numpy as np
import moose

chan = moose.HHChannel('channel')
chan.Xpower = 3
chan.Ypower = 1
m_gate = moose.element(f'{chan.path}/gateX')
h_gate = moose.element(f'{chan.path}/gateY')
vmin = -110e-3
vmax = 50e-3
vdivs = 1000
v = np.linspace(vmin, vmax, vdivs + 1) - (-70e-3)
m_alpha = (
    1e3 * (25 - v * 1e3) / (10 * (np.exp((25 - v * 1e3) / 10) - 1))
)
m_beta = 1e3 * 4 * np.exp(-v * 1e3 / 18)
m_gate.min = vmin
m_gate.max = vmax
m_gate.divs = vdivs
m_gate.tableA = m_alpha
m_gate.tableB = m_alpha + m_beta
h_alpha = 1e3 * 0.07 * np.exp(-v / 20e-3)
h_beta = 1e3 * 1 / (np.exp((30e-3 - v) / 10e-3) + 1)
h_gate.min = vmin
h_gate.max = vmax
h_gate.divs = vdivs
h_gate.tableA = h_alpha
h_gate.tableB = h_alpha + h_beta  
```

Q: How do I model an ion channel so that the channel equations are evaluated exactly instead of using lookup/interpolation tables?
----------------------------------------------------------------------------------------------------------------------------------

### A: Use the `HHChannelF` class for formula based evaluation of gating variables.

`HHChannelF` is similar to `HHChannel` but allows you to specify
expression strings for gate dynamics. It allows both alpha-beta form and
the tau-inf forms. To use alpha-beta form assign the formula for alpha
and beta to the fields `alpha` and `beta` of the gate in the channel.

Q: How do I model ion Hodgkin-Huxley-type channels that depend on both voltage and calcium concentration?
---------------------------------------------------------------------------------------------------------

### A: If your channel depends on two independent variables, for example voltage and calcium concentration, use `HHChannel2D` or `HHChannelF2D`.
