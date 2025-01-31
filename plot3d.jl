using PlotlyJS

# Generate some example data
x = 1:10
y = 1:10
z = [sin(xi) * cos(yi) for xi in x, yi in y]

# Create a 3D surface plot
trace = surface(x = x, y = y, z = z)
layout = Layout(title = "3D Surface Plot", scene = attr(xaxis_title = "X-axis", yaxis_title = "Y-axis", zaxis_title = "Z-axis"))
plot = Plot([trace], layout)

# Save the plot as an HTML file to preserve interactivity
savefig(plot, "3Dplot.html")

# Optionally, save the plot as a static image (e.g., PNG)
savefig(plot, "3Dplot.png")