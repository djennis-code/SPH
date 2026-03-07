using Raylib_cs;
using System.Numerics;

namespace SPH;

internal static class Program
{
    public static Vector2 windowSize = new(1200, 800);
    public static HydroDynamics simulation = new HydroDynamics();

    public static void InitSPH()
    {
        simulation.particleCount = 900;
        simulation.particleMass = 1.0f;
        simulation.particleRadius = 5.0f;
        simulation.bounds = windowSize;

        simulation.smoothingRadius = 75.0f;
        simulation.NormaliseKernels();

        simulation.gamma = 7f;
        simulation.refDensity = 1.2f;
        simulation.pressureConst = simulation.refDensity * 6f / simulation.gamma; // B = rho * cs^2 / gamma
        simulation.pressureExt = 1.0f;

        simulation.deltaTime = 0.05f;
        simulation.damping = 0.5f;
        simulation.viscosity = 0.5f;

        simulation.CreateBuffers();
        simulation.InitParticles();
    }

    public static void Main()
    {
        InitSPH();

        Raylib.InitWindow((int)windowSize.X, (int)windowSize.Y, "Smoothed Particle Magneto Dynamics");

        while (!Raylib.WindowShouldClose())
        {
            Raylib.BeginDrawing();

            Raylib.ClearBackground(Color.Black);

            simulation.ComputeDensities();
            simulation.ComputeAccelerations();
            simulation.Integrate();
            RenderParticles();

            Raylib.DrawFPS(5, 5);
            Raylib.EndDrawing();
        }

        Raylib.CloseWindow();
    }

    public static void RenderParticles()
    {
        int count = simulation.particleCount;
        float radius = simulation.particleRadius;
        Vector2[] positions = simulation.positions;

        for (uint i = 0; i < count; i++)
        {
            Color col = i != 550 ? Color.White : Color.Red;
            Vector2 pos = positions[i];
            Raylib.DrawCircle((int)pos.X, (int)pos.Y, radius, col);
        }
    }
}

