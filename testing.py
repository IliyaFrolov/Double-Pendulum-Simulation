class Vehicle():

    def __init__(self, wheels, colour):
        self.wheels = wheels
        self.colour = colour

    def change_colour(self, colour):
        self.colour = colour
    
    def print_random(self, something):
        print(something)

class Car(Vehicle):

    def __init__(self, wheels, colour):
        super().__init__(wheels, colour)
    
    def test(self, something):
        print(super())

toyota = Car(4, 'red')
toyota.test('hello')




