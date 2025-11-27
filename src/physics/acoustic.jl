
struct DirichletType <: FieldType end
#Dirchlet type of boundary conditions
struct NeumannType <: FieldType end
#Neumann type of boundary conditions


function greens(field::DirichletType, medium::Acoustic{2,T}, x::SVector{2,T}, ω::Float64) where T
   
    return (im/4) * hankelh1(0, ω / medium.c * norm(x))
end
