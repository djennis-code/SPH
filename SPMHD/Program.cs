using Raylib_cs;
using System.Numerics;

namespace SPH;

internal static class Program
{
    public static Vector2 windowSize = new(1280, 720);
    public static HydroDynamics simulation = new HydroDynamics();


    // Lambda = 1280 px, kappa = 2 Pi / Lambda ~ 0.005 / px.
    // Omega = Viscosity * kappa^2 ~ 2.5e-5 / px.

    public static void InitSPH()
    {
        simulation.particleCount = 900;
        simulation.particleMass = 1.0f;
        simulation.particleRadius = 3.0f;
        simulation.bounds = windowSize;
        simulation.wallForce = 100f;

        simulation.smoothingRadius = 50.0f;
        simulation.NormaliseKernels();

        simulation.gamma = 7f;
        simulation.refDensity = 0.3f;
        simulation.pressureConst = simulation.refDensity * 500.0f / simulation.gamma; // B = rho * cs^2 / gamma
        simulation.pressureExt = 0.0f;

        simulation.deltaTime = 0.03f; // seconds // frame
        simulation.viscosity = 1.0f; // Pa * frame

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

            simulation.IntegrateVerlet();
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
            Color col = Color.Blue;
            Vector2 pos = positions[i];
            Raylib.DrawCircle((int)pos.X, (int)pos.Y, radius, col);
        }
    }
}

