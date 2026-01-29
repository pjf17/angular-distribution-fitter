from math import sqrt

class Measurement:
    #class for handling error propagation
    def __init__(self, value, error):
        self.value = float(value)
        self.error = float(abs(error))  # Ensure error is positive
    
    def __str__(self):
        return f"{self.value:.3} ± {self.error:.3}"
    
    def __repr__(self):
        return f"Measurement({self.value}, {self.error})"
    
    def __add__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        new_value = self.value + other.value
        new_error = sqrt(self.error**2 + other.error**2)
        return Measurement(new_value, new_error)
    
    def __sub__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        new_value = self.value - other.value
        new_error = sqrt(self.error**2 + other.error**2)
        return Measurement(new_value, new_error)
    
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        new_value = self.value * other.value
        rel_error1 = self.error / abs(self.value) if self.value != 0 else 0
        rel_error2 = other.error / abs(other.value) if other.value != 0 else 0
        new_error = abs(new_value) * sqrt(rel_error1**2 + rel_error2**2)
        return Measurement(new_value, new_error)
    
    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        if other.value == 0:
            raise ValueError("Division by zero")
        new_value = self.value / other.value
        rel_error1 = self.error / abs(self.value) if self.value != 0 else 0
        rel_error2 = other.error / abs(other.value) if other.value != 0 else 0
        new_error = abs(new_value) * sqrt(rel_error1**2 + rel_error2**2)
        return Measurement(new_value, new_error)
    
    # To allow operations with numbers on the left (e.g., 2 * Measurement)
    def __radd__(self, other):
        return self.__add__(other)
    
    def __rsub__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        return other.__sub__(self)
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        return other.__truediv__(self)