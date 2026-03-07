using Raylib_cs;
using System.Numerics;

namespace SPH;

internal static class Program
{
    public static Vector2 windowSize = new(1200, 800);
    public static HydroDynamics simulation = new HydroDynamics();

    public static void InitSPH()
    {
        simulation.particleCount = 625;
        simulation.particleMass = 1.0f;
        simulation.particleRadius = 5.0f;
        simulation.bounds = windowSize;

        simulation.smoothingRadius = 80.0f;
        simulation.NormaliseKernels();

        simulation.refDensity = 1.5f;
        simulation.pressureConst = 1.5f; // B = rho * cs^2 / gamma
        simulation.pressureExt = 1.0f;

        simulation.deltaTime = 0.1f;
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

            simulation.DrawParticles();
            simulation.ComputeDensities();
            simulation.ComputeAccelerations();
            simulation.Integrate();

            //debug
            //int j = 50;
            //float rho = simulation.densities[j];
            //Vector2 vel = simulation.velocities[j];
            //Console.WriteLine(vel);

            Raylib.DrawFPS(5, 5);
            Raylib.EndDrawing();
        }

        Raylib.CloseWindow();
    }
}

