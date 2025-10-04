"""
Physical properties of the fluid
"""
@with_kw struct FluidProperties
    ρ::Float64          # Density [kg/m³]
    σ::Float64          # Surface tension [N/m]
    ν::Float64          # Kinematic viscosity [m²/s]
    ν_corrected::Float64 = 0.8025 * ν  # Corrected viscosity for Faraday
    g::Float64 = 9.8    # Gravity [m/s²]
end

"""
Physical properties of the solid
"""
@with_kw struct SolidProperties
    R₀::Float64         # Radius [m]
    ρ_s::Float64        # Density [kg/m³]
    m::Float64 = (4π/3) * R₀^3 * ρ_s  # Mass [kg]
end

"""
Dimensionless numbers based on characteristic scales L, V
"""
@with_kw struct DimensionlessNumbers
    Re::Float64         # Reynolds = LV/ν
    Fr::Float64         # Froude² = V²/(gL)
    We::Float64         # Weber = ρV²L/σ
    M::Float64          # Mass ratio = m/(ρL³)
    D::Float64 = 0.0    # Drag coefficient
end

"""
Compute dimensionless numbers for sphere impact
"""
function dimensionless_sphere_impact(fluid::FluidProperties, 
                                     solid::SolidProperties,
                                     V₀::Float64)
    L = solid.R₀
    V = V₀
    
    Re = L * V / fluid.ν
    Fr = V^2 / (fluid.g * L)
    We = fluid.ρ * V^2 * L / fluid.σ
    M = solid.m / (fluid.ρ * L^3)
    
    return DimensionlessNumbers(Re=Re, Fr=Fr, We=We, M=M, D=0.0)
end

"""
Compute dimensionless numbers for bouncing droplet (Faraday forcing)
"""
function dimensionless_bouncing_droplet(fluid::FluidProperties,
                                        solid::SolidProperties,
                                        f_forcing::Float64;
                                        use_corrected_viscosity::Bool=true)
    # Faraday wavelength (dispersion relation)
    ω = 2π * f_forcing
    # Solve: ω² = g*k + (σ/ρ)*k³ for k, then λ = 2π/k
    # Approximate for 20 cSt oil at 80 Hz:
    λ_F = 0.004969  # m (from paper)
    
    L = λ_F
    V = λ_F * f_forcing
    
    ν_use = use_corrected_viscosity ? fluid.ν_corrected : fluid.ν
    
    Re = L * V / ν_use
    Fr = (λ_F * f_forcing^2) / fluid.g
    We = fluid.ρ * λ_F^3 * f_forcing^2 / fluid.σ
    M = (4π * solid.R₀^3) / (3 * λ_F^3)
    
    # Stokes drag coefficient
    μ_air = 1.8e-5  # kg/(m·s)
    D = (9 * μ_air) / (2 * solid.ρ_s * solid.R₀^2 * f_forcing)
    
    return DimensionlessNumbers(Re=Re, Fr=Fr, We=We, M=M, D=D)
end

"""
Standard water properties at 25°C
"""
function water_properties()
    return FluidProperties(
        ρ = 1000.0,           # kg/m³
        σ = 0.0749,           # N/m
        ν = 8.94e-7           # m²/s
    )
end

"""
20 cSt silicone oil properties
"""
function silicone_oil_20cst()
    return FluidProperties(
        ρ = 949.0,            # kg/m³
        σ = 0.0206,           # N/m
        ν = 2.0e-5            # m²/s (20 cSt)
    )
end