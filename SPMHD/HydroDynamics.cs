using Raylib_cs;
using System;
using System.Formats.Asn1;
using System.Numerics;
using System.Runtime.InteropServices;

namespace SPH;

public struct Vector2i(int x, int y)
{
    public int X = x;
    public int Y = y;
}

public class HydroDynamics
{
    public const float PI = 3.141593f;
    
    public int particleCount = 25;
    public float particleMass = 1f;
    public float particleRadius = 1f;

    public Vector2 spawnOffset = new Vector2(0f, 0f);
    public Vector2 bounds = new(1f, 1f);
    public float domainThickness = 1f;
    public float gridSpacing = 1f;
    public float wallForce = 1f;

    public float smoothingRadius = 1e-2f;
    public float norm = 1f / (PI * 5000f);
    public float dNorm = 1f / (PI * 666.67f);

    public float gamma = 7f; // gamma = 7 for water
    public float refDensity = 1f;
    public float cSound = 1f;
    public float pressureConst = 1f;
    public float pressureExt = 1f;

    public float deltaTime = 1f;
    public float viscosity = 1f;
    public float gravity = 0f;

    public Vector2[] positions;
    public Vector2[] velocities;
    public Vector2[] accelerations;
    public float[] densities;

    public int gridLength = 1;
    public int gridHeight = 1;
    public List<int>[] grid;

    public Vector2i[] neighbours = new Vector2i[9];

    public void CreateBuffers()
    {
        positions = new Vector2[particleCount];
        velocities = new Vector2[particleCount];
        accelerations = new Vector2[particleCount];
        densities = new float[particleCount];
    }

    public void InitParticles()
    {
        int gridLength = (int) MathF.Sqrt(particleCount);

        for (uint i = 0; i < gridLength; i++)
        {
            for (uint j = 0; j < gridLength + 1; j++)
            {
                uint index = (uint) (i * gridLength + j);

                if (index >= particleCount) continue;

                Vector2 offset = new Vector2((bounds.X - gridSpacing * gridLength) / 2,
                                            (bounds.Y - gridSpacing * gridLength) / 2);

                positions[index] = new Vector2(gridSpacing * i + offset.X, gridSpacing * j + offset.Y) + spawnOffset;
            }
        }
    }

    public void InitHashMap()
    {
        gridLength = (int) (bounds.X / smoothingRadius) + 1;
        gridHeight = (int) (bounds.Y / smoothingRadius) + 1;

        int N = gridLength * gridHeight;
        grid = new List<int>[N];

        for (int n = 0; n < N; n++)
        {
            grid[n] = new List<int>();
        }

        uint i = 0;

        for (int x = -1; x <= 1; x++)
        {
            for (int y = -1; y <= 1; y++)
            {
                neighbours[i] = new Vector2i(x, y);
                i++;
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
        norm = 10f / (PI * MathF.Pow(smoothingRadius, 5));
        dNorm = 30f / (PI * MathF.Pow(smoothingRadius, 5));
    }

    public float DensityToPressure(float density)
    {
        float y = density / refDensity;
        return pressureConst * (MathF.Pow(y, gamma) - 1f) + pressureExt;
    }

    public float ComputeDensityAtPoint(Vector2 position)
    {
        float density = 0f;

        int cellX = (int) (position.X / smoothingRadius);
        int cellY = (int) (position.Y / smoothingRadius);

        foreach (Vector2i cellOffset in neighbours)
        {
            int nx = cellX + cellOffset.X;
            int ny = cellY + cellOffset.Y;

            if (nx < 0 || ny < 0 || nx >= gridLength || ny >= gridHeight)
                continue;

            int neighbour = ny * gridLength + nx;

            foreach (int j in grid[neighbour])
            {
                //if (i == j) continue;

                Vector2 offset = position - positions[j];
                float dist = offset.Length();
                float weight = Kernel(dist);

                density += weight * particleMass;
            }
        }
        return density;
    }

    public void ComputeDensities()
    {
        Parallel.For(0, particleCount, i =>
        {
            Vector2 pos = positions[i];
            densities[i] = ComputeDensityAtPoint(pos);
        });
    }

    public void ComputeAccelerations()
    {
        float s2 = smoothingRadius * smoothingRadius;
        float invTime = 1f / (deltaTime * deltaTime);

        Parallel.For(0, particleCount, i =>
        {
            Vector2 pos = positions[i];
            Vector2 vel = velocities[i];
            float rhoi = densities[i];

            Vector2 force = new();
            Vector2 laplacian = new(); // laplacian required for the viscosity contribution

            int cellX = (int)(pos.X / smoothingRadius);
            int cellY = (int)(pos.Y / smoothingRadius);

            foreach (Vector2i cellOffset in neighbours)
            {
                int nx = cellX + cellOffset.X;
                int ny = cellY + cellOffset.Y;

                if (nx < 0 || ny < 0 || nx >= gridLength || ny >= gridHeight)
                    continue;

                int neighbour = ny * gridLength + nx;

                foreach (int j in grid[neighbour])
                {
                    if (i == j) continue;

                    Vector2 offset = pos - positions[j];
                    Vector2 velOffset = vel - velocities[j];

                    float dist2 = offset.LengthSquared();
                    float rhoj = densities[j];

                    if (dist2 >= s2) continue;

                    float dist = MathF.Max(MathF.Sqrt(dist2), 0.01f);

                    Vector2 dir = offset / dist;
                    Vector2 gradient = DerivativeKernel(dist) * dir;

                    float sharedPressure = DensityToPressure(rhoi) / (rhoi * rhoi) + DensityToPressure(rhoj) / (rhoj * rhoj);

                    laplacian += 8f / rhoj * Vector2.Dot(offset, velOffset) / (dist2 + 0.01f * smoothingRadius) * gradient;
                    force += - sharedPressure * gradient;
                }
            }
            accelerations[i] = particleMass * (force + laplacian * viscosity) / rhoi;
        });
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
                acc.X -= wallForce * (smoothingRadius - dist.X) / particleMass;
            }
            else if (pos.X < smoothingRadius)
            {
                acc.X += wallForce * (smoothingRadius - pos.X) / particleMass;
            }

            if (dist.Y < smoothingRadius)
            {
                acc.Y -= wallForce * (smoothingRadius - dist.Y) / particleMass;
            }
            else if (pos.Y < smoothingRadius)
            {
                acc.Y += wallForce * (smoothingRadius - pos.Y) / particleMass;
            }

            accelerations[i] = acc;
        }
    }
    
    public void IntegrateEuler()
    {
        WriteHashMap();

        Vector2 g = new(0f, gravity);

        ComputeDensities();
        ComputeAccelerations();
        ApplyBoundaryConditions();

        for (uint i = 0; i < particleCount; i++)
        {   
            velocities[i] += (accelerations[i] + g) * deltaTime;
            positions[i] += velocities[i] * deltaTime;
        }

    }

    public void IntegrateVerlet()
    {
        Vector2 g = new(0f, gravity);

        for (uint i = 0; i < particleCount; i++)
        {
            positions[i] += velocities[i] * deltaTime + 0.5f * (accelerations[i] + g) * deltaTime * deltaTime;
        }

        WriteHashMap();

        Vector2[] prevAccelerations = (Vector2[]) accelerations.Clone();

        ComputeDensities();
        ComputeAccelerations();
        ApplyBoundaryConditions();

        for (uint i = 0; i < particleCount; i++)
        {
            velocities[i] += 0.5f * (prevAccelerations[i] + accelerations[i]) * deltaTime;
        }
    }

    public int GetCell(Vector2 pos)
    {
        int cellX = (int)(pos.X / smoothingRadius);
        int cellY = (int)(pos.Y / smoothingRadius);

        cellX = Math.Clamp(cellX, 0, gridLength - 1);
        cellY = Math.Clamp(cellY, 0, gridHeight - 1);

        return cellY * gridLength + cellX;
    }

    public void WriteHashMap()
    {
        // Clear grid
        for (int i = 0; i < grid.Length; i++)
        {
            grid[i].Clear();
        }

        // Insert particles
        for (int i = 0; i < particleCount; i++)
        {
            int cell = GetCell(positions[i]);
            grid[cell].Add(i);
        }
    }
}


