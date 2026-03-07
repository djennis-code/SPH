using Raylib_cs;
using System;
using System.Formats.Asn1;
using System.Numerics;
using System.Runtime.InteropServices;

namespace SPH;

public class HydroDynamics
{
    public const float PI = 3.141593f;
    
    public int particleCount = 25;
    public float particleMass = 1f;
    public float particleRadius = 1f;

    public Vector2 bounds = new(100f, 100f);
    public int boundaryParticles = 0;
    public float wallForce = 10f;

    public float smoothingRadius = 10f;
    public float norm = 1f / (PI * 5000f);
    public float dNorm = 1f / (PI * 666.67f);

    public float gamma = 7f; // gamma = 7 for water
    public float refDensity = 1f;
    public float pressureConst = 1f;
    public float pressureExt = 1f;

    public float deltaTime = 1f;
    public float damping = 1f;
    public float viscosity = 1f;

    public Vector2[] positions;
    public Vector2[] velocities;
    public Vector2[] accelerations;
    public float[] densities;

    public Vector2[] boundaryPositions;

    public void CreateBuffers()
    {
        positions = new Vector2[particleCount];
        velocities = new Vector2[particleCount];
        accelerations = new Vector2[particleCount];
        densities = new float[particleCount];
    }

    public void InitBoundary()
    {
        int nx = (int) (bounds.X / smoothingRadius);
        int ny = (int) (bounds.Y / smoothingRadius);

        boundaryParticles = 2 * nx + 2 * ny;
        boundaryPositions = new Vector2[boundaryParticles];

        for (int x = 0; x < nx; x++)
        {
            Vector2 pos1 = new(x * smoothingRadius, 0f);
            Vector2 pos2 = new(x * smoothingRadius, bounds.Y);

            boundaryPositions[x + 0] = pos1;
            boundaryPositions[x + 1] = pos2;
        }

        for (int y = 0; y < ny; y++)
        {
            Vector2 pos1 = new(0f, y * smoothingRadius);
            Vector2 pos2 = new(bounds.X, y * smoothingRadius);

            boundaryPositions[2 * nx + y + 0] = pos1;
            boundaryPositions[2 * nx + y + 1] = pos2;
        }
    }

    public void InitParticles()
    {
        int gridLength = (int) MathF.Sqrt(particleCount);
        float gridSpacing = 5f * particleRadius;

        for (uint i = 0; i <= gridLength; i++)
        {
            for (uint j = 0; j <= gridLength + 1; j++)
            {
                uint index = (uint)(i * gridLength + j);

                if (index >= particleCount) continue;

                Vector2 offset = new Vector2((bounds.X - gridSpacing * gridLength) / 2,
                                            (bounds.Y - gridSpacing * gridLength) / 2);

                positions[index] = new Vector2(gridSpacing * i + offset.X, gridSpacing * j + offset.Y);
            }
        }
    }

    public float Kernel(float dist)
    {
        float x = Math.Max(smoothingRadius - dist, 0f);
        return x * x * x * norm;
    }

    public float DerivativeKernel(float dist)
    {
        float x = Math.Max(smoothingRadius - dist, 0f);
        return - x * x * dNorm;
    }

    public void NormaliseKernels()
    {
        norm = 15f / (PI * MathF.Pow(smoothingRadius, 4));
        dNorm = 45f / (PI * MathF.Pow(smoothingRadius, 3));
    }

    public float DensityToPressure(float density)
    {
        float y = density / refDensity;
        return pressureConst * (MathF.Pow(y, gamma) - 1f) + pressureExt;
    }

    public void ComputeDensities()
    {
        for (uint i = 0; i < particleCount; i++)
        {
            Vector2 pos = positions[i];
            float density = 0f;

            for (uint j = 0; j < particleCount; j++)
            {
                //if (i == j) continue;

                Vector2 offset = pos - positions[j];
                float dist = offset.Length();
                float weight = Kernel(dist);

                density += weight * particleMass;   
            }

            densities[i] = density;
        }
    }

    public void ComputeAccelerations()
    {
        Vector2 g = new Vector2(0f, 2.0f);
        float s2 = smoothingRadius * smoothingRadius;

        for (uint i = 0; i < particleCount; i++)
        {
            Vector2 pos = positions[i];
            Vector2 vel = velocities[i];
            float rhoi = densities[i];

            Vector2 force = new();
            Vector2 vis = new(); // force vector for the viscosity contribution

            for (uint j = 0; j < particleCount; j++)
            {
                if (i == j) continue;

                Vector2 offset = pos - positions[j];
                float dist2 = offset.LengthSquared();
                float rhoj = densities[j];

                if (dist2 >= s2) continue;

                float dist = MathF.Max(MathF.Sqrt(dist2), 0.1f);
                Vector2 dir = offset / dist;

                float weight = Kernel(dist);
                float gradient = DerivativeKernel(dist);

                float sharedPressure = (DensityToPressure(rhoi) / (rhoi * rhoi)  + DensityToPressure(rhoj) / (rhoj * rhoj));

                vis += weight * particleMass * (vel - velocities[j]) / rhoj;
                force += - dir * gradient * sharedPressure * particleMass;
            }   
            accelerations[i] = (force - vis * viscosity + g) / rhoi;
        }
    }

    public void ApplyBoundaryConditions()
    {
        for (uint i = 0; i < particleCount; i++)
        {
            Vector2 pos = positions[i];
            Vector2 dist = bounds - pos;
            Vector2 acc = accelerations[i];

            if (dist.X < smoothingRadius)
            {
                acc.X -= wallForce * (smoothingRadius - dist.X) * (smoothingRadius - dist.X);
            }
            if (pos.X < smoothingRadius)
            {
                acc.X += wallForce * (smoothingRadius - pos.X) * (smoothingRadius - pos.X);
            }

            if (dist.Y < smoothingRadius)
            {
                acc.Y -= wallForce * (smoothingRadius - dist.Y) * (smoothingRadius - dist.Y);
            }
            if (pos.Y < smoothingRadius)
            {
                acc.Y += wallForce * (smoothingRadius - pos.Y) * (smoothingRadius - pos.Y);
            }

            accelerations[i] = acc;
        }
    }
    
    public void Integrate()
    {
        ApplyBoundaryConditions();

        for (uint i = 0; i < particleCount; i++)
        {
            // Eulerian integration
            
            velocities[i] += accelerations[i] * deltaTime;
            positions[i] += velocities[i] * deltaTime;
        }

    }
}


