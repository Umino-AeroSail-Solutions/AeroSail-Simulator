import pygame
import numpy as np
import Sail as S
import Profile as P

P.initializeXfoil('C:/Xfoil699src', 'C:/Xfoil699src/xfoil.exe')
sail_instance = S.Sail_Class('Data/E473coordinates.txt', 5, 0.4, 30, panels=20)
sail_instance.load_interpolation('Data/interpolationCR4sail_XFLR5.npz')
interpolation = 'Data/interpolationCR4sail_XFLR5.npz'

# Initialize Pygame
pygame.init()

# Screen dimensions
WIDTH, HEIGHT = 800, 600
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Sail and Flap Visualizer")

# Font setup
pygame.font.init()
font = pygame.font.SysFont('Arial', 24)

# Colors
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
RED = (255, 0, 0)
BLUE = (0, 0, 255)
GREEN = (0, 255, 0)
YELLOW = (255, 255, 0)

# Sail parameters
sail_length = 120
flap_length = 80
thrust_length = 100
sail_angle = 0
flap_angle = 0


def draw_sail(screen, origin, angle, color):
    end_pos = (origin[0] + sail_length * np.cos(angle), origin[1] - sail_length * np.sin(angle))
    pygame.draw.line(screen, color, origin, end_pos, 5)
    return end_pos


def draw_flap(screen, origin, angle, color):
    end_pos = (origin[0] + flap_length * np.cos(angle), origin[1] - flap_length * np.sin(angle))
    pygame.draw.line(screen, color, origin, end_pos, 3)


def draw_thrust(screen, origin, magnitude, color):
    # Thrust vector is horizontal (right direction)
    end_pos = (origin[0] - magnitude, origin[1])
    pygame.draw.line(screen, color, origin, end_pos, 3)
    pygame.draw.polygon(screen, color, [
        (end_pos[0], end_pos[1]),
        (end_pos[0] + 10, end_pos[1] - 5),
        (end_pos[0] + 10, end_pos[1] + 5)
    ])


def draw_awa_arrow(screen, origin, angle, length, color):
    end_pos = (origin[0] + length * np.cos(angle), origin[1] - length * np.sin(angle))
    pygame.draw.line(screen, color, origin, end_pos, 3)
    pygame.draw.polygon(screen, color, [
        (end_pos[0], end_pos[1]),
        (end_pos[0] - 10 * np.cos(angle - np.pi / 6), end_pos[1] + 10 * np.sin(angle - np.pi / 6)),
        (end_pos[0] - 10 * np.cos(angle + np.pi / 6), end_pos[1] + 10 * np.sin(angle + np.pi / 6))
    ])

def draw_arrow(screen, origin, angle, length, color):
    end_pos = (origin[0] + length * np.cos(angle), origin[1] - length * np.sin(angle))
    pygame.draw.line(screen, color, origin, end_pos, 3)
    pygame.draw.polygon(screen, color, [
        (end_pos[0], end_pos[1]),
        (end_pos[0] - 10 * np.cos(angle - np.pi / 6), end_pos[1] + 10 * np.sin(angle - np.pi / 6)),
        (end_pos[0] - 10 * np.cos(angle + np.pi / 6), end_pos[1] + 10 * np.sin(angle + np.pi / 6))
    ])


def draw_text(screen, text, position, color):
    text_surface = font.render(text, True, color)
    screen.blit(text_surface, position)


# Main loop
running = True
while running:
    screen.fill(BLACK)
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    mx, my = pygame.mouse.get_pos()
    AWA = (np.arctan2((HEIGHT // 2) - my, mx - (WIDTH // 2))) + np.pi
    sail_instance.get_opt_pos(AWA)
    sail_angle = AWA - np.radians(sail_instance.opt_alpha)
    flap_angle = sail_angle - sail_instance.opt_flap  # Example flap adjustment, change as needed

    # Calculate thrust vector magnitude cf
    cf = sail_instance.cts.max()  # Replace with actual calculation logic
    # Draw sail, flap, thrust vector, and AWA arrow
    sail_origin = (WIDTH // 2, HEIGHT // 2)
    sail_end = draw_sail(screen, sail_origin, sail_angle, BLUE)
    draw_flap(screen, sail_end, flap_angle, RED)
    draw_thrust(screen, sail_origin, cf * 100, GREEN)
    draw_awa_arrow(screen, (mx,my), AWA, 100, YELLOW)

    # Draw text for angles in degrees and thrust magnitude
    sail_angle_degrees = sail_instance.opt_alpha
    flap_angle_degrees = np.degrees(sail_instance.opt_flap)
    coefficients = sail_instance.get_sail_coefficients(abs(sail_instance.opt_alpha), abs(sail_instance.opt_flap), s_interpolation=interpolation)
    cl = coefficients[0]
    cd = coefficients[1]

    draw_text(screen, f"Sail Angle: {sail_angle_degrees:.2f}°", (20, 20), WHITE)
    draw_text(screen, f"Flap Angle: {flap_angle_degrees:.2f}°", (20, 50), WHITE)
    draw_text(screen, f"Thrust Coefficient Magnitude: {cf:.2f}", (20, 80), WHITE)
    draw_text(screen, f"Lift Coefficient Magnitude: {cl:.2f}", (20, 110), WHITE)
    draw_text(screen, f"Drag Coefficient Magnitude: {cd:.2f}", (20, 140), WHITE)
    draw_text(screen, f"Cl/Cd: {(cl/cd):.2f}", (20, 170), WHITE)

    # Draw lift
    if sail_angle_degrees > 0:
        draw_arrow(screen, sail_origin, AWA+(np.pi/2), cl*100, (0,10,100))
    else:
        draw_arrow(screen, sail_origin, AWA-(np.pi/2), cl*100, (0,10,100))
    # Draw drag
    draw_arrow(screen, sail_origin, AWA, cd * 100, (100, 10, 0))

    # Update the display
    pygame.display.flip()

# Quit Pygame
pygame.quit()
