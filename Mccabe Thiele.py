import numpy as np
import matplotlib.pyplot as plt

# Parameters
def main():
    print("\nInput parameters and press Enter to continue")

    try:
        alpha = float(input("Alpha (relative volatility): "))
        molar_reflux_ratio = float(input("Molar Reflux Ratio: "))
        feed_thermal_condition = float(input("Feed Thermal Condition (q) (1=Saturated Liquid): "))
        zf = float(input("Feed Composition (zf): "))
        xb = float(input("Bottoms Composition (xb): "))
        xd = float(input("Distillate Composition (xd): "))
    except ValueError:
        print("Error: All inputs must be valid numbers.")
        return main()  # restart if conversion fails

    # Basic validity checks
    if alpha <= 1:
        print("Error: Relative volatility (alpha) must be greater than 1.")
        return main()

    if molar_reflux_ratio < 0:
        print("Error: Molar reflux ratio must be non-negative.")
        return main()

    if not (0 <= feed_thermal_condition <= 2):
        print("Error: Feed thermal condition (q) must be between 0 and 2.")
        return main()

    if not (0 <= xb < zf < xd <= 1):
        print("Error: Must satisfy 0 ≤ xb < zf < xd ≤ 1 for valid compositions.")
        return main()

    print("\n Inputs received successfully:")
    print(f"Alpha = {alpha}")
    print(f"Molar Reflux Ratio = {molar_reflux_ratio}")
    print(f"Feed Thermal Condition (q) = {feed_thermal_condition}")
    print(f"Feed Composition (zf) = {zf}")
    print(f"Bottoms Composition (xb) = {xb}")
    print(f"Distillate Composition (xd) = {xd}")

    # Grids
    x = np.linspace(0, 1, 500)

    # Equilibrium relationship
    def y_eq(x):
        return (alpha * x) / (1 + (alpha - 1) * x)

    def x_from_y(y):
        denom = (alpha - y * (alpha - 1))
        return y / denom

    # Rectifying line
    m_rect = molar_reflux_ratio / (molar_reflux_ratio + 1)
    b_rect = 1 / (molar_reflux_ratio + 1)
    def y_rect(x):
        return m_rect * x + b_rect

    # q-line
    if abs(feed_thermal_condition - 1.0) < 1e-6:
        q_vertical = True
        x_qline = zf
    else:
        q_vertical = False
        q_slope = feed_thermal_condition / (feed_thermal_condition - 1)
        q_intercept = -zf / (feed_thermal_condition - 1)
        def y_qline(x):
            return q_slope * x + q_intercept

    # Intersection of rectifying line and q-line
    if q_vertical:
        x_intersect = x_qline
        y_intersect = y_rect(x_intersect)
    else:
        x_intersect = (q_intercept - b_rect) / (m_rect - q_slope)
        y_intersect = y_rect(x_intersect)

    # Stripping line (corrected)
    m_strip = (y_intersect - xb) / (x_intersect - xb)
    b_strip = y_intersect - m_strip * x_intersect
    def y_strip(x):
        return m_strip * x + b_strip

    x_vals, y_vals = [], []
    stages = 0
    y_current = xd

    for _ in range(200):
        # Horizontal step (to eq curve)
        x_eq_point = x_from_y(y_current)
        x_vals += [x_eq_point, x_eq_point]
        y_vals += [y_current, y_current]
        stages += 0.5

        # Vertical step (to operating line)
        if q_vertical:
            if x_eq_point >= x_qline:
                y_next = y_rect(x_eq_point)
            else:
                y_next = y_strip(x_eq_point)
        else:
            if x_eq_point >= x_intersect:
                y_next = y_rect(x_eq_point)
            else:
                y_next = y_strip(x_eq_point)

        x_vals[-1] = x_eq_point
        y_vals[-1] = y_next
        stages += 0.5
        y_current = y_next

        if x_eq_point <= xb + 1e-3:
            break

    n_stages = int(np.ceil(stages))
    # --- Plotting ---
    plt.figure(figsize=(8, 8))

    # Diagonal and equilibrium curve
    plt.plot(x, x, label='y = x (Diagonal)', color='blue')
    plt.plot(x, y_eq(x), label='Equilibrium Curve', color='orange')

    # Operating lines
    plt.plot(x, y_rect(x), label='Rectifying Line', color='green')
    plt.plot(x, y_strip(x), label='Stripping Line', color='red')
    if q_vertical:
        plt.axvline(x=x_qline, linestyle="--", color='gray', label="q-line (vertical)")
    else:
        plt.plot(x, y_qline(x), linestyle="--", color='gray', label="q-line")

    # Stages
    plt.step(x_vals, y_vals, where='post', label='Stages', color='purple')

    # Markers for key points
    plt.scatter([zf], [zf], color='red', marker='x', s=80, label='Feed (zf)')
    plt.scatter([xd], [xd], color='green', marker='o', s=80, label='Distillate (xd)')
    plt.scatter([xb], [xb], color='blue', marker='s', s=80, label='Bottoms (xb)')

    # Labels and formatting
    plt.xlabel('x (Liquid Mole Fraction)')
    plt.ylabel('y (Vapor Mole Fraction)')
    plt.title('McCabe-Thiele Diagram')
    plt.legend()
    plt.grid(True)
    plt.axis("square")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    print("Save as image if needed press s.")
    plt.show()
    if input().lower() == 's': 
        plt.savefig("mccabe_thiele_diagram.png")
    else:
        plt.close()
    

    
    print(f"Number of stages required: {n_stages}")


if __name__ == "__main__":
    main()
    