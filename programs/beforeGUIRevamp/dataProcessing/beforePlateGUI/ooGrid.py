root = Tk()

fill = 1

class MyButton(Button): #Convention is for class names to start with uppercase letters
    def __init__(self, master):
        super(MyButton, self).__init__(master, image = images[0], borderwidth = 0)
        self.type = 0
        self.already_changed = False

    def change(self):
        if self.type == fill:
            self.type = 0
        else:
            self.type = fill
        self.config(image=images[self.type])

    def mouse_entered(self):
        if not self.already_changed:
            self.change()
            self.already_changed = True

    def mouse_up(self):
        self.already_changed = False

class Container(Frame):
    def __init__(self, master, width, height):
        super(Container, self).__init__(master)

        buttons = []

        for y in range(height):
            buttons.append([])
            for x in range(width):
                button = MyButton(self)
                button.grid(row = x, column = y)

                buttons[y].append(button)

        self.buttons = buttons

        self.bind_all("<Button-1>", self.mouse_down)
        self.bind_all("<ButtonRelease-1>", self.mouse_up)
        self.bind_all("<B1-Motion>", self.mouse_motion)

        self.mouse_pressed = False

    def mouse_down(self, e):
        self.update_containing_button(e)
        self.mouse_pressed = True

    def mouse_up(self, e):
        self.mouse_pressed = False
        for row in self.buttons:
            for button in row:
                button.mouse_up()

    def mouse_motion(self, e):
        self.update_containing_button(e)

    def update_containing_button(self, e):
        for row in self.buttons:
            for button in row:
                if self.winfo_containing(e.x_root, e.y_root) is button:
                    button.mouse_entered()

grid = Container(root, 15, 15)
grid.pack()

root.mainloop()
