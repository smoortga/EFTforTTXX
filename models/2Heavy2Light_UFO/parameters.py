# This file was automatically created by FeynRules 2.3.28
# Mathematica version: 10.2.0 for Mac OS X x86 (64-bit) (July 7, 2015)
# Date: Tue 7 Nov 2017 10:51:54



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
cabi = Parameter(name = 'cabi',
                 nature = 'external',
                 type = 'real',
                 value = 0.227736,
                 texname = '\\theta _c',
                 lhablock = 'CKMBLOCK',
                 lhacode = [ 1 ])

C83Qq = Parameter(name = 'C83Qq',
                  nature = 'external',
                  type = 'real',
                  value = 1,
                  texname = '\\frac{\\text{Subsuperscript}\\left[C,q_L Q_L,83\\right]}{\\Lambda ^2}',
                  lhablock = 'DIM6',
                  lhacode = [ 1 ])

C81Qq = Parameter(name = 'C81Qq',
                  nature = 'external',
                  type = 'real',
                  value = 1,
                  texname = '\\frac{\\text{Subsuperscript}\\left[C,q_L Q_L,81\\right]}{\\Lambda ^2}',
                  lhablock = 'DIM6',
                  lhacode = [ 2 ])

C13Qq = Parameter(name = 'C13Qq',
                  nature = 'external',
                  type = 'real',
                  value = 1,
                  texname = '\\frac{\\text{Subsuperscript}\\left[C,q_L Q_L,13\\right]}{\\Lambda ^2}',
                  lhablock = 'DIM6',
                  lhacode = [ 3 ])

C11Qq = Parameter(name = 'C11Qq',
                  nature = 'external',
                  type = 'real',
                  value = 1,
                  texname = '\\frac{\\text{Subsuperscript}\\left[C,q_L Q_L,11\\right]}{\\Lambda ^2}',
                  lhablock = 'DIM6',
                  lhacode = [ 4 ])

C8td = Parameter(name = 'C8td',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\frac{\\text{Subsuperscript}\\left[C,d_R t_R,1\\right]}{\\Lambda ^2}',
                 lhablock = 'DIM6',
                 lhacode = [ 5 ])

C8tu = Parameter(name = 'C8tu',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\frac{\\text{Subsuperscript}\\left[C,t_R u_R,1\\right]}{\\Lambda ^2}',
                 lhablock = 'DIM6',
                 lhacode = [ 6 ])

C1td = Parameter(name = 'C1td',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\frac{\\text{Subsuperscript}\\left[C,d_R t_R,1\\right]}{\\Lambda ^2}',
                 lhablock = 'DIM6',
                 lhacode = [ 7 ])

C1tu = Parameter(name = 'C1tu',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\frac{\\text{Subsuperscript}\\left[C,t_R u_R,1\\right]}{\\Lambda ^2}',
                 lhablock = 'DIM6',
                 lhacode = [ 8 ])

C8Qd = Parameter(name = 'C8Qd',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\frac{\\text{Subsuperscript}\\left[C,d_R Q_L,8\\right]}{\\Lambda ^2}',
                 lhablock = 'DIM6',
                 lhacode = [ 9 ])

C8Qu = Parameter(name = 'C8Qu',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\frac{\\text{Subsuperscript}\\left[C,Q_L u_R,8\\right]}{\\Lambda ^2}',
                 lhablock = 'DIM6',
                 lhacode = [ 10 ])

C8tq = Parameter(name = 'C8tq',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\frac{\\text{Subsuperscript}\\left[C,q_L t_R,8\\right]}{\\Lambda ^2}',
                 lhablock = 'DIM6',
                 lhacode = [ 11 ])

C1Qd = Parameter(name = 'C1Qd',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\frac{\\text{Subsuperscript}\\left[C,d_R Q_L,1\\right]}{\\Lambda ^2}',
                 lhablock = 'DIM6',
                 lhacode = [ 12 ])

C1Qu = Parameter(name = 'C1Qu',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\frac{\\text{Subsuperscript}\\left[C,Q_L u_R,1\\right]}{\\Lambda ^2}',
                 lhablock = 'DIM6',
                 lhacode = [ 13 ])

C1tq = Parameter(name = 'C1tq',
                 nature = 'external',
                 type = 'real',
                 value = 1,
                 texname = '\\frac{\\text{Subsuperscript}\\left[C,q_L t_R,1\\right]}{\\Lambda ^2}',
                 lhablock = 'DIM6',
                 lhacode = [ 14 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.0000116637,
               texname = 'G_f',
               lhablock = 'SMINPUTS',
               lhacode = [ 2 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.1184,
               texname = '\\alpha _s',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

ymb = Parameter(name = 'ymb',
                nature = 'external',
                type = 'real',
                value = 4.7,
                texname = '\\text{ymb}',
                lhablock = 'YUKAWA',
                lhacode = [ 5 ])

ymt = Parameter(name = 'ymt',
                nature = 'external',
                type = 'real',
                value = 172,
                texname = '\\text{ymt}',
                lhablock = 'YUKAWA',
                lhacode = [ 6 ])

ymtau = Parameter(name = 'ymtau',
                  nature = 'external',
                  type = 'real',
                  value = 1.777,
                  texname = '\\text{ymtau}',
                  lhablock = 'YUKAWA',
                  lhacode = [ 15 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 172,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.7,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

MH = Parameter(name = 'MH',
               nature = 'external',
               type = 'real',
               value = 125,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.4952,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.085,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.50833649,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

WH = Parameter(name = 'WH',
               nature = 'external',
               type = 'real',
               value = 0.00407,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 25 ])

CKM1x1 = Parameter(name = 'CKM1x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.cos(cabi)',
                   texname = '\\text{CKM1x1}')

CKM1x2 = Parameter(name = 'CKM1x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.sin(cabi)',
                   texname = '\\text{CKM1x2}')

CKM2x1 = Parameter(name = 'CKM2x1',
                   nature = 'internal',
                   type = 'complex',
                   value = '-cmath.sin(cabi)',
                   texname = '\\text{CKM2x1}')

CKM2x2 = Parameter(name = 'CKM2x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.cos(cabi)',
                   texname = '\\text{CKM2x2}')

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\alpha _{\\text{EW}}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

MW = Parameter(name = 'MW',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2))))',
               texname = 'M_W')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')

sw2 = Parameter(name = 'sw2',
                nature = 'internal',
                type = 'real',
                value = '1 - MW**2/MZ**2',
                texname = '\\text{sw2}')

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - sw2)',
               texname = 'c_w')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(sw2)',
               texname = 's_w')

g1 = Parameter(name = 'g1',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g_1')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

vev = Parameter(name = 'vev',
                nature = 'internal',
                type = 'real',
                value = '(2*MW*sw)/ee',
                texname = '\\text{vev}')

lam = Parameter(name = 'lam',
                nature = 'internal',
                type = 'real',
                value = 'MH**2/(2.*vev**2)',
                texname = '\\text{lam}')

yb = Parameter(name = 'yb',
               nature = 'internal',
               type = 'real',
               value = '(ymb*cmath.sqrt(2))/vev',
               texname = '\\text{yb}')

yt = Parameter(name = 'yt',
               nature = 'internal',
               type = 'real',
               value = '(ymt*cmath.sqrt(2))/vev',
               texname = '\\text{yt}')

ytau = Parameter(name = 'ytau',
                 nature = 'internal',
                 type = 'real',
                 value = '(ymtau*cmath.sqrt(2))/vev',
                 texname = '\\text{ytau}')

muH = Parameter(name = 'muH',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(lam*vev**2)',
                texname = '\\mu')

I1a33 = Parameter(name = 'I1a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb',
                  texname = '\\text{I1a33}')

I2a33 = Parameter(name = 'I2a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I2a33}')

I3a33 = Parameter(name = 'I3a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I3a33}')

I4a33 = Parameter(name = 'I4a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb',
                  texname = '\\text{I4a33}')

