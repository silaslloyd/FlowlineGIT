function IfThenElse(condition,t,f) 
   if condition then
     return t
   else
     return f
   end 
end

function smoothplus(x, eps)
  k = 10
  X = eps+ math.log(1+math.exp(k*x))/k
  return X
end

function smoothplusDerivative(x)
  k = 10
  return 1/(1+math.exp(k*x))
end


function smoothmax(a,b)
  k = 50
  return 1/k * math.log(math.exp(k*a) + math.exp(k*b))
end

function smoothmaxDerivative(a,b)
  k = 50
  return math.exp(k*a)/(math.exp(k*a) + math.exp(k*b))
end

function hwCondition(gm)
  if gm < 0.5 then 
    return hr
  else 
    return 0
  end
end

function Tinit(x,z,Tsurf,Tdivide,T0,xg,N_total,N_Bedrock)

  zmax = 1
  zmin = 0.5
  z = (z-zmin)/(zmax-zmin)
  Tb = Tbed(x,Tdivide,xg)
  T = Tb + (Tsurf - Tb)*z
  T = IfThenElse(T > 273.15, 273.15, T)
  return T
end

function Tbed(x,Tdivide,xg)
    Tb = Tdivide + x/(xg) *(T0-Tdivide)
    Tb = IfThenElse(Tb > 273.15, 273.15, Tb)
    return Tb
end

function Tdeep(x,z,Tsurf,Tdivide,T0,xg,N_total,N_Bedrock)
  zmax= 0.5
  zmin = 0
  z =((z-zmin)/(zmax-zmin)) -1
  Tbase = 273.15 + 10000*gthf/heatconductivityrock
  Tb =  Tbed(x,Tdivide,xg)
  T = Tb - z*(Tbase - Tb)
  return T
end


function H(x)
  return IfThenElse((x > 0.0),6.0*(math.sqrt(x + 5.0e3) - math.sqrt(5.0e3)),0.0) + Hmin
end

function mboscillations(t)
  period = 1000
  return aDIM*(1+cos(2*3.14*t/period))/2
end

function coupleflux(gwflx,x)
  if (x > 0.0) then
    return yearinsec*gwflx
  else
    return 0.0
  end  
end

function waterpressure(p,x)
  if (x > 0) then
    return p*MPainPa
  else
    return p0
  end
end  

function watersource(x)
  relx = x - 50000.0
  dist = math.abs(relx)
  deltax = 30000.0
  if (dist < deltax) then
    return 0.9*(1.0 - dist/deltax)
  else  
    return 0.0
  end  
end

function waterpressure_SI(p,x)
  if (x > 0) then
    return p
  else
    return p0
  end
end

function TDepSlidingExp(C,Th,gamma)
  Th = Th - 273.15
  Th = IfThenElse(Th > 0.0, 0.0, Th)
  factor = 1.0 --(mask+1.0)/2 --IfThenElse(mask<=0.0,0.000001,1.0)
  return factor*C / math.exp(Th/gamma)
end

function TDepSlidingCoulomb(C, N, As, q, n, xvelo, yvelo, zvelo, gamma, Th, normalstress)
  Th = Th - 273.15
  Th = IfThenElse(Th > 0.0, 0.0, Th)
  -- basal velocity 
  -- since bed is impermeable the basal velocity is the velocity magnitude at the bed 
  ub = math.sqrt(xvelo*xvelo + yvelo*yvelo + zvelo*zvelo)
  a = (q-1.0)^(q-1.0)/(q^q)
  chi = ub/(C^n * math.min(N, normalstress)^n *As) -- make sure that possible negative water pressure doesnt play a role
  Cw = 1/((As*(1 + a*chi^q))^(1.0/n))
  return Cw / math.exp(Th/gamma)
end

function getHorizontalAbsoluteVelocity(vx, vy, mask)
  if (mask > 0) then
    return math.sqrt(vx*vx + vy+vy)
  else
    return 0.0
  end  
end


function getAbsoluteVelocity(vx, vy, vz)
    return math.sqrt(vx*vx + vy*vy + vz*vz)
end

function bedrock1(x)
  return (720.0 - (720.0 + 372.54)*x/1.0526e+06) + 200.0
end

function bedrock(x)
  return 720.0 - 778.5*x/750.0e3 + 0.1
end  

function settimestepsize(nt)
  dt = dtInit*(dtIncr^nt)
  print(">>>>>>>> Current timestep",dt)
  if ( dt > dtMax ) then
    dt = dtMax
  end
  return dt
end

function min(a,b )
 if (a > b) then
   return b
 else
   return a
 end  
end

function max(a,b )
 if (a < b) then
   return b
 else
   return a
 end  
end

function InitialSheet(x)
  return 1e-4 + 0*0.05*math.exp(-((x-3e6)/5e4)^2)
end

function getmeltrate(loads,weights,friction, gm, time)
  if gm <= 0.0 then
    return 0.0
  else
    melt =  -(loads-0*friction)/(weights*rhoi*Lw)
    melt = IfThenElse(time < 1E-2, 0, melt)
    --melt = -(loads-friction)/(weights)
    return melt
  end
end

function InitialCTMask(T,T0)
  
  return IfThenElse(T <272, -1, 1)
end

function getmeltrate2(lx, ly, lz, vx, vy, vz)
 frictionheating = -1.0*(min(lx*vx,0.0) + min(ly*vy,0.0) + min(lz*vz,0.0))
 -- print(">>>>>>>> Friction Heating:",frictionheating/(rhoi*Lw))
 return frictionheating/(rhoi*Lw)
 -- return 1.0/(365.25*24*3600) 
 -- return meltrate
end

function arrhenius(Trel)
  if (Trel > -10.0) then
    A0 = A2
    Q = Q2
  else
   A0 = A1
    Q = Q1
  end
  T = min(0.0,Trel)
  A = A0 * math.exp(-Q/(Rg*(273.15 - T)))
  -- print(">>>>>>>> Arrhenius", T, A)
  return A 
end

function sheetthickness(N, vx, vy, vz, T)
  ub = math.max(math.sqrt(vx*vx + vy*vy + vz*vz), 1e-14)
  factor = lr*N/(ub*viscosity) + 1
  hw = hr/factor
  --if (T == 0.0) then
  --  hw = hr/factor 
  --else
  --  hw = 1e-5
  --end
  return hw
end 

function sheetconductivity(N, vx, vy, vz, potentialgradientabs, T)
  hw = math.max(sheetthickness(N, vx, vy, vz, T), 1e-14)--sheetthickness(N, vx, vy, vz)--math.max(sheetthickness(N, vx, vy, vz), 1e-25)
  gradient = math.max(potentialgradientabs, 1e-14)
  conductivity = k*hw^(1.25)/(math.sqrt(gradient))
  return conductivity
end

function frictionloads(lx,ly,lz,vx,vy,vz)
  cx = max(-lx*vx,0)
  cy = max(-ly*vy,0)
  cz = max(-lz*vz,0)
  return cx + cy + cz
end  

function max(a,b)
  if (a > b) then
    return a
  else
    return b
  end
end