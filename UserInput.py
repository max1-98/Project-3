import pygame
import sys

# Initialize Pygame
pygame.init()

# Set up the display
width, height = 400, 300
screen = pygame.display.set_mode((width, height))
pygame.display.set_caption("Pygame Buttons")

# Define colors
white = (255, 255, 255)
black = (0, 0, 0)
gray = (200, 200, 200)

# Define fonts
font = pygame.font.Font(None, 36)

# Define button properties
button_width, button_height = 150, 50
button_padding = 20

# Define button positions
button1_pos = ((width - button_width) // 2, (height - button_height) // 2 - button_padding - button_height)
button2_pos = ((width - button_width) // 2, (height - button_height) // 2 + button_padding)

# Define button rectangles
button1_rect = pygame.Rect(button1_pos, (button_width, button_height))
button2_rect = pygame.Rect(button2_pos, (button_width, button_height))

# Define states
MAIN_MENU = 0
LINEAR_ODE_INPUT = 1
NON_LINEAR_ODE_INPUT = 2

# Initial state
current_state = MAIN_MENU


# Define input box properties
input_box_width, input_box_height = 200, 40
input_box_pos = ((width - input_box_width) // 2, height // 2 - input_box_height // 2)

# Define button properties
button_width, button_height = 150, 50
button_pos = ((width - button_width) // 2, height // 2 + input_box_height)

# Define button rectangle
button_rect = pygame.Rect(button_pos, (button_width, button_height))

# Input variables
input_text = ''
input_active = False

def ninput():

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()
        elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:  # Left mouse button
            if button_rect.collidepoint(event.pos):
                try:
                    n = int(input_text)
                    print(f"The order of the ODE is {n}")
                except ValueError:
                    print("Please enter a valid integer.")
        elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:  # Left mouse button
            if button_rect.collidepoint(event.pos):
                input_active = not input_active
            else:
                input_active = False
        elif event.type == pygame.KEYDOWN:
            if input_active:
                if event.key == pygame.K_RETURN:
                    try:
                        n = int(input_text)
                        print(f"The order of the ODE is {n}")
                    except ValueError:
                        print("Please enter a valid integer.")
                elif event.key == pygame.K_BACKSPACE:
                    input_text = input_text[:-1]
                else:
                    input_text += event.unicode

    # Draw the background
    screen.fill(white)

    # Draw the input box
    pygame.draw.rect(screen, gray, (input_box_pos[0], input_box_pos[1], input_box_width, input_box_height))
    input_surface = font.render(input_text, True, black)
    screen.blit(input_surface, (input_box_pos[0] + 5, input_box_pos[1] + 5))

    # Draw the button
    pygame.draw.rect(screen, gray, button_rect)
    button_text = font.render("Submit", True, black)
    screen.blit(button_text, (button_rect.centerx - button_text.get_width() // 2, button_rect.centery - button_text.get_height() // 2))

    # Update the display
    pygame.display.flip()

    # Control the frame rate
    pygame.time.Clock().tick(60)








# Run the game loop
while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()
        elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:  # Left mouse button
            if button1_rect.collidepoint(event.pos):
                current_state = LINEAR_ODE_INPUT
            elif button2_rect.collidepoint(event.pos):
                current_state = NON_LINEAR_ODE_INPUT

    # Draw the background
    screen.fill(white)

    # Draw the buttons
    pygame.draw.rect(screen, gray, button1_rect)
    pygame.draw.rect(screen, gray, button2_rect)

    # Draw the text on the buttons
    button1_text = font.render("Linear ODE", True, black)
    button2_text = font.render("Non-Linear ODE", True, black)
    screen.blit(button1_text, (button1_rect.centerx - button1_text.get_width() // 2, button1_rect.centery - button1_text.get_height() // 2))
    screen.blit(button2_text, (button2_rect.centerx - button2_text.get_width() // 2, button2_rect.centery - button2_text.get_height() // 2))

    if current_state == LINEAR_ODE_INPUT:
        ninput()
        pass
    elif current_state == NON_LINEAR_ODE_INPUT:
        # Code for Non-Linear ODE input page
        pass
    

    # Update the display
    pygame.display.flip()

    # Control the frame rate
    pygame.time.Clock().tick(60)
