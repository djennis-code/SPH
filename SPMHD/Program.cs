using Raylib_cs;
using System.Numerics;

namespace SPH;

internal static class Program
{
    public static Vector2 domainSize = new(2.0f, 1.0f); //  meter

    public static float upScale = 800f; //  pixels / meter
    public static Vector2 windowSize = upScale * domainSize;

    public static HydroDynamics simulation = new();

    public static void InitSPH()
    {
        simulation.particleCount = 36*36;
        simulation.smoothingRadius = 8e-2f; //  meter
        simulation.NormaliseKernels();

        simulation.gamma = 7f; //  exponent
        simulation.refDensity = 1e3f; //  Kg / meter^2
        simulation.cSound = 1.5e1f; //  meter / frame
        simulation.pressureConst = simulation.refDensity * simulation.cSound * simulation.cSound / simulation.gamma; // Pascal; B = rho * cs^2 / gamma
        simulation.pressureExt = 1e5f; //  Pascal

        simulation.particleMass = simulation.refDensity * (simulation.smoothingRadius * simulation.smoothingRadius * 3.14f); //  Kg; m = rho * V
        simulation.particleRadius = 4.0f; //  pixels
        simulation.bounds = domainSize;

        simulation.wallForce = 1e8f; //  meter / frame^2
        simulation.gridSpacing = 0.2f * simulation.smoothingRadius; //  meter
        simulation.gravity = 1e4f; //  meter / frame^2

        simulation.deltaTime = 1e-4f; //  seconds / frame
        simulation.viscosity = 5e6f; //  Pascal * frame

        simulation.CreateBuffers();
        simulation.InitParticles();
        simulation.InitHashMap();
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
            Vector2 pos = upScale * positions[i];

            Raylib.DrawCircle((int)pos.X, (int)pos.Y, radius, col);
        }
    }
}

